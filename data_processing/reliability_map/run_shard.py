#!/usr/bin/env python3
"""
run_shard.py — score one shard of the pooled SGE variants with AlphaGenome.

You have several API keys. Run this once per key, in parallel, each with a
different --shard-id, so each key is responsible for ~1/num_shards of the
35,333 variants. Each process writes its own per-variant wide table; combine
them afterwards with combine_shards.py.

Example (6 keys, 6 shards, run these 6 lines in 6 terminals / processes):
    python run_shard.py --api-key $KEY0 --num-shards 6 --shard-id 0
    python run_shard.py --api-key $KEY1 --num-shards 6 --shard-id 1
    ...
    python run_shard.py --api-key $KEY5 --num-shards 6 --shard-id 5

Resumable: if interrupted, just rerun the same command — variants already in the
shard's output file are skipped. One progress line is printed per variant.

The API key is read only from --api-key (or $ALPHAGENOME_API_KEY) and is never
written to any output file.
"""
import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ag_common as ag  # noqa: E402

from alphagenome.models import dna_client  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--api-key", default=os.environ.get("ALPHAGENOME_API_KEY"),
                   help="AlphaGenome API key (or set $ALPHAGENOME_API_KEY). "
                        "Never hard-code this into the script.")
    p.add_argument("--num-shards", type=int, required=True,
                   help="Total number of shards (= number of API keys you run).")
    p.add_argument("--shard-id", type=int, required=True,
                   help="This process's shard index, 0..num_shards-1.")
    p.add_argument("--input", default=str(HERE.parent / "analysis" / "pooled_labeled.csv"),
                   help="Path to pooled_labeled.csv.")
    p.add_argument("--out-dir", default=str(HERE / "out"),
                   help="Directory for per-shard output tables.")
    p.add_argument("--seq-length", default=ag.DEFAULT_SEQ_LENGTH,
                   choices=list(ag.SEQUENCE_LENGTHS.keys()),
                   help="Sequence context window centred on the variant "
                        "(smaller = faster/cheaper; default 1MB = full context).")
    p.add_argument("--max-workers", type=int, default=4,
                   help="Concurrent in-flight requests for this one key.")
    p.add_argument("--checkpoint-every", type=int, default=50,
                   help="Rewrite the shard output every N completed variants.")
    return p.parse_args()


def atomic_write_csv(df, path):
    tmp = str(path) + ".tmp"
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def main():
    args = parse_args()
    if not args.api_key:
        sys.exit("ERROR: provide --api-key or set $ALPHAGENOME_API_KEY")
    if not (0 <= args.shard_id < args.num_shards):
        sys.exit(f"ERROR: --shard-id must be in 0..{args.num_shards - 1}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"shard_{args.shard_id:02d}_of_{args.num_shards:02d}.csv"

    # Strided shard selection keeps every gene represented in every shard.
    df = ag.load_variants(args.input)
    shard = df.iloc[args.shard_id::args.num_shards].reset_index(drop=True)

    # Resume: skip variants already scored in this shard's output.
    done = set()
    rows = []
    if out_path.exists():
        prev = pd.read_csv(out_path)
        rows = prev.to_dict("records")
        done = set(prev["accession"].astype(str))
    todo = shard[~shard["accession"].astype(str).isin(done)].reset_index(drop=True)

    print(f"[shard {args.shard_id}/{args.num_shards}] {len(shard)} variants total, "
          f"{len(done)} already done, {len(todo)} to score, seq_length={args.seq_length}, "
          f"workers={args.max_workers}", flush=True)
    if len(todo) == 0:
        print(f"[shard {args.shard_id}/{args.num_shards}] nothing to do -> {out_path}", flush=True)
        return

    client = dna_client.create(args.api_key)
    scorers = ag.get_scorers()

    def work(rec):
        variant = ag.make_variant(rec)
        tidy = ag.score_one(client, variant, scorers, seq_length_key=args.seq_length)
        return ag.reduce_to_wide(tidy, rec["gene"], str(rec["accession"]))

    t0 = time.time()
    n_ok = n_err = 0
    records = todo.to_dict("records")
    with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
        futs = {ex.submit(work, r): r for r in records}
        for i, fut in enumerate(as_completed(futs), 1):
            rec = futs[fut]
            acc = str(rec["accession"])
            try:
                rows.append(fut.result())
                n_ok += 1
                status = rows[-1].get("ag_status", "ok")
            except Exception as e:  # noqa: BLE001
                rows.append({"accession": acc, "ag_gene": rec["gene"],
                             "ag_status": f"error: {type(e).__name__}: {e}"})
                n_err += 1
                status = "ERROR"
            rate = i / max(time.time() - t0, 1e-9)
            print(f"[shard {args.shard_id}/{args.num_shards}] "
                  f"{i}/{len(records)} {acc} ({rec['gene']}) {status} "
                  f"[{rate:.2f} var/s]", flush=True)
            if i % args.checkpoint_every == 0:
                atomic_write_csv(pd.DataFrame(rows), out_path)

    atomic_write_csv(pd.DataFrame(rows), out_path)
    dt = time.time() - t0
    print(f"[shard {args.shard_id}/{args.num_shards}] DONE {n_ok} ok, {n_err} errors, "
          f"{dt/60:.1f} min -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
