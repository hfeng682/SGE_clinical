#!/usr/bin/env python3
"""
build_candidates.py -- rebuild the candidate set from the pooled screen table.

The candidates are the variants where DNA-based practice would see a protein-coding
change (missense) or no change (synonymous), but our screens' mRNA score shows the
transcript drops (rna_drop5 == True). These are the variants whose splicing
consequence is hidden behind a protein-coding annotation -- the natural place to
look for splice effects that sequence-based tools miss.

Reads : ../data/pooled_labeled.csv   (the frozen pooled screen table)
Writes: data/candidate_variants.csv

Offline. Run once; run_spliceai.py and analyze.py consume the output.
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
POOLED = os.path.join(HERE, "..", "data", "pooled_labeled.csv")
OUT = os.path.join(HERE, "data", "candidate_variants.csv")

KEEP = ["gene", "chrom", "pos", "ref", "alt", "hgvs_nt", "consequence",
        "coarse_consequence", "route_class", "protein_visible",
        "dist_to_junction", "junction_side", "dist_bin",
        "score", "call", "rna_score", "rna_measured", "rna_drop5", "rna_drop1"]


def main():
    df = pd.read_csv(POOLED)
    cand = df[df.coarse_consequence.isin(["missense_variant", "synonymous_variant"])
              & (df.rna_drop5 == True)].copy()
    cand = cand[[c for c in KEEP if c in cand.columns]].copy()
    cand["chrom"] = cand.chrom.astype(str)
    cand["variant_hg38"] = cand.apply(
        lambda r: f"{'chr'+str(r.chrom) if not str(r.chrom).startswith('chr') else r.chrom}"
                  f"-{r.pos}-{r.ref}-{r.alt}", axis=1)
    cand = cand.sort_values(["gene", "pos"]).reset_index(drop=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    cand.to_csv(OUT, index=False)
    print(f"[build_candidates] {len(cand)} candidates -> {OUT}")
    print("  by gene:", cand.gene.value_counts().to_dict())
    print("  by consequence:", cand.coarse_consequence.value_counts().to_dict())
    print("  by distance bin:", cand.dist_bin.value_counts().to_dict())


if __name__ == "__main__":
    main()
