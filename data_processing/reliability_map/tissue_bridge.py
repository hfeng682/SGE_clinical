#!/usr/bin/env python3
"""
tissue_bridge.py — AlphaGenome tissue-transfer evidence for UNMEASURABLE variants.

The SGE assay cannot measure mRNA abundance for intronic variants. Where the
reliability map proved AG reliable (the junction-proximal zone), AG supplies tissue-
resolved evidence for those variants. Everything deeper is flagged unreliable
and given no AG evidence.

Reliability gate (strict, from the reliability map): only the junction-proximal distance
bins (analysis_common.RELIABLE_DIST_BINS) are gated in. This mirrors the one
reliable reliability map stratum (exonic edge <=3 nt) by distance symmetry.

Output (results/):
  tissue_bridge.csv   one row per unmeasurable variant: reliability flag, AG splice
                      call, disease-tissue + pan-tissue expression evidence,
                      per-GTEx-tissue vector, and a transfer verdict.

Run:  python tissue_bridge.py
"""
import os
import sys
import json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analysis_common as C


def disease_evidence(row, gtex_cols, short):
    tissues = C.DISEASE_TISSUES.get(row.gene, [])
    cols = [c for c in gtex_cols if short[c] in tissues]
    vals = row[cols].astype(float).values
    vals = vals[~np.isnan(vals)]
    if len(vals) == 0:
        return pd.Series({"disease_expr_lfc_mean": np.nan,
                          "disease_expr_lfc_maxabs": np.nan, "n_disease_tissues": 0})
    return pd.Series({"disease_expr_lfc_mean": vals.mean(),
                      "disease_expr_lfc_maxabs": vals[np.argmax(np.abs(vals))],
                      "n_disease_tissues": len(vals)})


def main():
    os.makedirs(C.RESULTS_DIR, exist_ok=True)
    df = C.load_table()
    u = C.unmeasurable(df)
    print(f"unmeasurable variants: {len(u)}")

    gtex_cols = sorted([c for c in u.columns if c.startswith("expr_lfc_gtex__")])
    short = {c: c.replace("expr_lfc_gtex__", "") for c in gtex_cols}

    # reliability gate (strict distance)
    u["ag_reliable"] = u[C.DIST_BIN].isin(C.RELIABLE_DIST_BINS)

    # AG splice call at the reliability map 90%-specificity threshold (recompute on measurable negatives)
    m = C.measurable(df).copy()
    m["ag_splice_mag"] = C.ag_splice_mag(m)
    neg = m.loc[~m[C.EVENT].astype(bool), "ag_splice_mag"].values
    thr90 = C.threshold_at_specificity(neg, 0.90)
    u["ag_splice_mag"] = C.ag_splice_mag(u)
    u["ag_splice_disruptor"] = u["ag_splice_mag"] >= thr90

    # expression evidence
    de = u.apply(lambda r: disease_evidence(r, gtex_cols, short), axis=1)
    u = pd.concat([u, de], axis=1)
    allv = u[gtex_cols].astype(float)
    u["pantissue_expr_lfc_mean"] = allv.mean(axis=1)
    u["tissue_lfc_std"] = allv.std(axis=1)   # cross-tissue spread: tissue-specificity of the effect

    # transfer verdict
    def verdict(r):
        if not r.ag_reliable:
            return "AG_unreliable_no_evidence"
        if not r.ag_splice_disruptor:
            return "no_splice_disruption_predicted"
        return ("splice_disruptor_expr_concordant" if r.disease_expr_lfc_mean < 0
                else "splice_disruptor_expr_discordant")
    u["transfer_verdict"] = u.apply(verdict, axis=1)

    front = ["accession", "gene", "chrom", "pos", "ref", "alt", "consequence",
             "coarse_consequence", "route_class", C.DIST, C.DIST_BIN, "junction_side",
             "ag_reliable", "ag_splice_mag", "ag_splice_disruptor",
             "splice_sites_maxabs", "splice_usage_maxabs", "splice_junctions_maxabs",
             C.EXPR_COL, C.EXPR_SIGNED_COL,
             "disease_expr_lfc_mean", "disease_expr_lfc_maxabs", "n_disease_tissues",
             "pantissue_expr_lfc_mean", "tissue_lfc_std", "transfer_verdict"]
    out = u[front + gtex_cols].copy()
    out.to_csv(os.path.join(C.RESULTS_DIR, "tissue_bridge.csv"), index=False)

    # console report
    gated = int(u.ag_reliable.sum())
    print(f"threshold (90% spec) = {thr90:.4f}")
    print(f"AG-reliable (gated in): {gated} / {len(u)} ({gated/len(u):.1%})")
    print("gated-in by dist_bin:")
    print(u[u.ag_reliable][C.DIST_BIN].value_counts().to_string())
    print("\ntransfer_verdict (gated-in):")
    print(u[u.ag_reliable].transfer_verdict.value_counts().to_string())
    disr = u[u.ag_reliable & u.ag_splice_disruptor]
    if len(disr):
        print(f"\nsplice-disruptors: {len(disr)}; disease-tissue expr negative "
              f"{(disr.disease_expr_lfc_mean < 0).mean():.0%}; "
              f"median cross-tissue std {disr.tissue_lfc_std.median():.4f} (near tissue-invariant)")


if __name__ == "__main__":
    main()
