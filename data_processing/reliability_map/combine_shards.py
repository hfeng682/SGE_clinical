#!/usr/bin/env python3
"""
combine_shards.py — merge the per-shard AlphaGenome tables into the single
end-product per-variant wide table, joined back to pooled_labeled.csv.

Run after all shards from run_shard.py have finished:
    python combine_shards.py --num-shards 6

It reads ./out/shard_XX_of_YY.csv, concatenates them, checks that every variant
in pooled_labeled.csv is covered exactly once, joins the AlphaGenome columns
onto the full labelled table, and writes:

    ./alphagenome_per_variant.csv   — one row per variant: all pooled_labeled.csv
                                      columns + AlphaGenome splice summaries and
                                      tissue-resolved expression columns.

Coverage and any scoring errors are reported so nothing is silently dropped.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--num-shards", type=int, required=True)
    p.add_argument("--input", default=str(HERE.parent / "analysis" / "pooled_labeled.csv"))
    p.add_argument("--shard-dir", default=str(HERE / "out"))
    p.add_argument("--out", default=str(HERE / "alphagenome_per_variant.csv"))
    p.add_argument("--allow-partial", action="store_true",
                   help="Write the table even if some variants are missing/errored "
                        "(default: refuse so an incomplete run is caught).")
    return p.parse_args()


def main():
    args = parse_args()
    shard_dir = Path(args.shard_dir)

    frames, found = [], []
    for sid in range(args.num_shards):
        path = shard_dir / f"shard_{sid:02d}_of_{args.num_shards:02d}.csv"
        if path.exists():
            frames.append(pd.read_csv(path))
            found.append(sid)
        else:
            print(f"WARNING: missing shard file {path}", flush=True)
    if not frames:
        sys.exit(f"ERROR: no shard files found in {shard_dir}")
    print(f"[combine] found shards {found} ({len(frames)}/{args.num_shards})", flush=True)

    ag = pd.concat(frames, ignore_index=True)
    dups = ag["accession"].duplicated().sum()
    if dups:
        print(f"WARNING: {dups} duplicate accessions across shards; keeping first "
              f"non-error occurrence.", flush=True)
        ag = ag.sort_values("ag_status").drop_duplicates("accession", keep="first")

    status = ag["ag_status"].astype(str)
    n_ok = int((status == "ok").sum())
    n_empty = int((status == "empty").sum())
    n_err = int(status.str.startswith("error").sum())
    print(f"[combine] AlphaGenome rows: {len(ag)}  (ok={n_ok}, empty={n_empty}, error={n_err})", flush=True)
    if n_err:
        errs = ag.loc[status.str.startswith("error"), ["accession", "ag_status"]].head(10)
        print("[combine] first errored variants:\n" + errs.to_string(index=False), flush=True)

    labelled = pd.read_csv(args.input, low_memory=False)
    covered = set(ag["accession"].astype(str))
    all_acc = set(labelled["accession"].astype(str))
    missing = all_acc - covered
    extra = covered - all_acc
    print(f"[combine] coverage: {len(all_acc & covered)}/{len(all_acc)} variants scored; "
          f"{len(missing)} missing, {len(extra)} extra (not in input)", flush=True)

    if (missing or n_err) and not args.allow_partial:
        sys.exit(f"ERROR: {len(missing)} variants missing and {n_err} errored. "
                 f"Rerun the relevant shards (resumable), or pass --allow-partial "
                 f"to write anyway.")

    # Left-join AlphaGenome columns onto the full labelled table (input is authoritative
    # for coverage). Drop the redundant helper column from the shard side.
    ag_join = ag.drop(columns=[c for c in ["ag_gene"] if c in ag.columns])
    merged = labelled.merge(ag_join, on="accession", how="left")

    # Order: keep input columns first, then AlphaGenome columns.
    ag_cols = [c for c in merged.columns if c not in labelled.columns]
    merged = merged[list(labelled.columns) + ag_cols]

    merged.to_csv(args.out, index=False)
    print(f"[combine] wrote {len(merged)} rows x {merged.shape[1]} cols -> {args.out}", flush=True)
    print(f"[combine] AlphaGenome columns added ({len(ag_cols)}): "
          f"{', '.join(ag_cols[:12])}{' ...' if len(ag_cols) > 12 else ''}", flush=True)


if __name__ == "__main__":
    main()
