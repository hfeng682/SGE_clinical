#!/usr/bin/env python3
"""
pilot_run.py — score ~50 variants with AlphaGenome so we can sanity-check the
output BEFORE committing to the full 35,333-variant run.

Run this once on your machine with any one of your API keys:
    python pilot_run.py --api-key $YOUR_KEY

It selects a stratified ~50-variant sample spanning all variant classes and both
score behaviours (the motivating missense-LoF-with-mRNA-drop cohort, protein-
blind synonymous/splice-region, the nonsense positive control, intronic, UTR),
scores them, and writes two files into ./pilot_out/:

    pilot_wide.csv  — the reduced per-variant wide table (what the full run emits)
    pilot_tidy.csv  — AlphaGenome's full tidy long output (one row per track),
                      so we can inspect the raw scores and tissue tracks

Send both files back and I'll check that AlphaGenome's splice/expression scores
track the measured RNA readout before you launch the sharded full run.

The API key is read only from --api-key (or $ALPHAGENOME_API_KEY) and is never
written to any output file.
"""
import argparse
import os
import sys
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ag_common as ag  # noqa: E402

from alphagenome.models import dna_client, variant_scorers  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--api-key", default=os.environ.get("ALPHAGENOME_API_KEY"),
                   help="AlphaGenome API key (or set $ALPHAGENOME_API_KEY).")
    p.add_argument("--input", default=str(HERE.parent / "analysis" / "pooled_labeled.csv"))
    p.add_argument("--out-dir", default=str(HERE / "pilot_out"))
    p.add_argument("--seq-length", default=ag.DEFAULT_SEQ_LENGTH,
                   choices=list(ag.SEQUENCE_LENGTHS.keys()))
    p.add_argument("--n-per-group", type=int, default=6,
                   help="Target variants per mechanism group (~48 total).")
    p.add_argument("--seed", type=int, default=7)
    return p.parse_args()


def select_pilot(df, n, seed):
    """Stratified ~50-variant sample spanning the mechanism contrasts."""
    def take(sub):
        return sub.sample(min(n, len(sub)), random_state=seed) if len(sub) else sub

    groups = [
        ("missense_LoF_rna_drop",
         df[(df.coarse_consequence == "missense_variant") & (df.call == "LoF") & (df.rna_drop5 == True)]),
        ("missense_LoF_no_drop",
         df[(df.coarse_consequence == "missense_variant") & (df.call == "LoF") & (df.rna_drop5 == False) & (df.rna_measured == True)]),
        ("missense_Normal",
         df[(df.coarse_consequence == "missense_variant") & (df.call == "Normal") & (df.rna_measured == True)]),
        ("synonymous",
         df[(df.coarse_consequence == "synonymous_variant") & (df.rna_measured == True)]),
        ("splice_region",
         df[(df.coarse_consequence.isin(["splicing_variant", "splice_site_variant"])) & (df.rna_measured == True)]),
        ("nonsense",
         df[df.coarse_consequence == "stop_gained"]),
        ("intron",
         df[df.coarse_consequence == "intron_variant"]),
        ("UTR",
         df[df.coarse_consequence == "UTR_variant"]),
    ]
    pieces = []
    for label, sub in groups:
        s = take(sub).copy()
        s["pilot_group"] = label
        pieces.append(s)
    return pd.concat(pieces).drop_duplicates("accession").reset_index(drop=True)


def main():
    args = parse_args()
    if not args.api_key:
        sys.exit("ERROR: provide --api-key or set $ALPHAGENOME_API_KEY")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = ag.load_variants(args.input)
    # columns select_pilot needs beyond the key cols
    for c in ["coarse_consequence", "call", "rna_measured", "rna_drop5", "rna_score", "dist_bin"]:
        if c not in df.columns:
            sys.exit(f"ERROR: input missing column needed for pilot selection: {c}")
    pilot = select_pilot(df, args.n_per_group, args.seed)
    print(f"[pilot] selected {len(pilot)} variants across "
          f"{pilot.pilot_group.nunique()} groups, {pilot.gene.nunique()} genes, "
          f"seq_length={args.seq_length}", flush=True)

    client = dna_client.create(args.api_key)
    scorers = ag.get_scorers()

    wide_rows, tidy_frames = [], []
    t0 = time.time()
    records = pilot.to_dict("records")
    for i, rec in enumerate(records, 1):
        acc = str(rec["accession"])
        variant = ag.make_variant(rec)
        try:
            tidy = ag.score_one(client, variant, scorers, seq_length_key=args.seq_length)
            wide = ag.reduce_to_wide(tidy, rec["gene"], acc)
            tidy = tidy.copy()
            tidy.insert(0, "accession", acc)
            tidy_frames.append(tidy)
            status = wide.get("ag_status", "ok")
        except Exception as e:  # noqa: BLE001
            wide = {"accession": acc, "ag_gene": rec["gene"],
                    "ag_status": f"error: {type(e).__name__}: {e}"}
            status = "ERROR"
        # carry the measured-readout columns for side-by-side inspection
        for c in ["gene", "coarse_consequence", "call", "rna_score", "rna_measured",
                  "rna_drop5", "dist_bin", "pilot_group"]:
            wide[f"measured_{c}"] = rec.get(c)
        wide_rows.append(wide)
        rate = i / max(time.time() - t0, 1e-9)
        print(f"[pilot] {i}/{len(records)} {acc} ({rec['gene']}, {rec['coarse_consequence']}) "
              f"{status} [{rate:.2f} var/s]", flush=True)

    wide_df = pd.DataFrame(wide_rows)
    wide_df.to_csv(out_dir / "pilot_wide.csv", index=False)
    if tidy_frames:
        pd.concat(tidy_frames, ignore_index=True).to_csv(out_dir / "pilot_tidy.csv", index=False)
    dt = time.time() - t0
    print(f"[pilot] DONE {len(wide_rows)} variants in {dt/60:.1f} min -> "
          f"{out_dir/'pilot_wide.csv'} and pilot_tidy.csv", flush=True)


if __name__ == "__main__":
    main()
