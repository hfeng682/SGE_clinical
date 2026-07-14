#!/usr/bin/env python3
"""
Reporting and quality-control for the clinical control set. Reads the files
produced by build_clinical_controls.py and reproduces every number quoted in
STEP0_REPORT.md, plus the summary figure.

Run after the build, from inside the controls/ directory:

    python step0_analysis.py

It prints:
    - the ClinVar/assay join rates (overall, SNV, indel) and the indel reconciliation;
    - the assay-vs-truth concordance table (sensitivity / specificity);
    - preliminary (unbounded) per-gene OddsPath and the implied tier;
    - the pooled per-stratum control counts.

And writes:
    control_summary.png
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils import paths
DERIV = paths.CLINVAR_GNOMAD
GENES = ["BARD1", "PALB2", "RAD51D", "BRCA1", "VHL"]

# OddsPath strength thresholds, pathogenic direction (Brnich et al. 2020).
ODDS_SUPPORTING, ODDS_MODERATE, ODDS_STRONG = 2.1, 4.3, 18.7


def tier(op: float) -> str:
    if not np.isfinite(op):          return "n/a"
    if op >= ODDS_STRONG:            return "Strong"
    if op >= ODDS_MODERATE:          return "Moderate"
    if op >= ODDS_SUPPORTING:        return "Supporting"
    return "below"


def load():
    ctrl = pd.read_csv(paths.CLINICAL_CONTROLS)
    ctrl["chrom"] = ctrl.chrom.astype(str)
    return ctrl


# join rates + indel reconciliation
def join_diagnostics(ctrl: pd.DataFrame) -> None:
    m = ctrl.clin_label.notna()
    is_indel = (ctrl.ref.str.len() != 1) | (ctrl.alt.str.len() != 1)
    print("ClinVar/assay join")
    print(f"  overall match : {int(m.sum()):>6} / {len(ctrl)}  ({100*m.mean():.1f}%)")
    print(f"  SNV match     : {100*m[~is_indel].mean():.1f}%")
    print(f"  indel match   : {100*m[is_indel].mean():.1f}%  "
          f"({int(m[is_indel].sum())} / {int(is_indel.sum())})")

    # 2-star pathogenic SNVs that DID join, as a fraction of those present in-window
    cv = pd.read_csv(DERIV / "clinvar_labeled_regions.csv").fillna({"gold_stars": 0})
    cv["chrom"] = cv.chrom.astype(str)
    cv_indel = (cv.ref.astype(str).str.len() != 1) | (cv.alt.astype(str).str.len() != 1)
    p2_snv = cv[(cv.gold_stars >= 2) & (cv.clin_label == "P/LP") & (~cv_indel)]
    akey = set(zip(ctrl.chrom, ctrl.pos, ctrl.ref, ctrl.alt))
    joined = sum((c, p, r, a) in akey for c, p, r, a in
                 zip(p2_snv.chrom, p2_snv.pos, p2_snv.ref, p2_snv.alt))
    print(f"  2-star P/LP SNVs joined: {joined} / {len(p2_snv)} "
          f"({100*joined/len(p2_snv):.1f}%) — shortfall is screen coverage, not coords")


# concordance: assay call vs clinical truth (2-star controls)
def concordance(ctrl: pd.DataFrame) -> pd.DataFrame:
    c = ctrl[(ctrl.gold_stars >= 2) & ctrl.clin_label.isin(["P/LP", "B/LB"])]
    ct = pd.crosstab(c.clin_label, c.call).reindex(
        index=["P/LP", "B/LB"], columns=["LoF", "Uncertain", "Normal"]).fillna(0).astype(int)
    P, B = c[c.clin_label == "P/LP"], c[c.clin_label == "B/LB"]
    print("\nConcordance (2-star controls)")
    print(ct.to_string())
    print(f"  LoF-call sensitivity : {(P.call=='LoF').mean():.3f}  "
          f"({int((P.call=='LoF').sum())}/{len(P)})")
    print(f"  Normal-call specificity: {(B.call=='Normal').mean():.3f}  "
          f"({int((B.call=='Normal').sum())}/{len(B)})")
    return ct


# preliminary (unbounded) OddsPath, pathogenic direction, per gene
def prelim_oddspath(ctrl: pd.DataFrame) -> None:
    c = ctrl[(ctrl.gold_stars >= 2) & ctrl.clin_label.isin(["P/LP", "B/LB"])]
    print("\nPreliminary (unbounded) OddsPath (PS3 direction), per gene")
    for g, sub in c.groupby("gene"):
        P1 = (sub.clin_label == "P/LP").mean()               # prior among controls
        lof = sub[sub.call == "LoF"]
        if len(lof) == 0 or P1 in (0, 1):
            print(f"  {g:7} n/a"); continue
        P2 = (lof.clin_label == "P/LP").mean()               # posterior among LoF-called
        if P2 >= 1:                                          # cap to avoid divide-by-zero
            P2 = 1 - 1 / (2 * len(lof))
        op = (P2 * (1 - P1)) / ((1 - P2) * P1)
        print(f"  {g:7} n_ctrl={len(sub):4d}  P1={P1:.2f}  OddsPath={op:6.1f}  -> {tier(op)}")
    print("  (point estimates on the full control set; no held-out split, no bound)")


# summary figure
def figure(ctrl: pd.DataFrame, ct: pd.DataFrame) -> None:
    c = ctrl[(ctrl.gold_stars >= 2) & ctrl.clin_label.isin(["P/LP", "B/LB"])]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))

    pg = pd.read_csv(DERIV / "overlap_per_gene.csv")
    pg = pg[pg.gene != "POOLED"]
    x, w = np.arange(len(pg)), 0.38
    ax[0].bar(x - w/2, pg.P_2star, w, label="P/LP (2★+)", color="#c0392b")
    ax[0].bar(x + w/2, pg.B_2star, w, label="B/LB (2★+)", color="#2e86c1")
    ax[0].set_xticks(x); ax[0].set_xticklabels(pg.gene); ax[0].set_ylabel("controls (n)")
    ax[0].set_title("a  Per-gene ClinVar controls (2★+)"); ax[0].legend(frameon=False)

    ctn = ct.div(ct.sum(1), axis=0)
    ax[1].imshow(ctn.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax[1].set_xticks(range(3)); ax[1].set_xticklabels(["LoF", "Uncertain", "Normal"])
    ax[1].set_yticks(range(2)); ax[1].set_yticklabels(["P/LP", "B/LB"])
    for i in range(2):
        for j in range(3):
            ax[1].text(j, i, f"{ct.values[i,j]}\n{ctn.values[i,j]*100:.0f}%",
                       ha="center", va="center",
                       color="white" if ctn.values[i, j] > 0.5 else "black")
    ax[1].set_title("b  Assay call vs ClinVar truth"); ax[1].set_xlabel("assay call")

    st = pd.read_csv(DERIV / "overlap_per_stratum.csv")
    pool = st.groupby("route_class").agg(P2=("P2", "sum"), B2=("B2", "sum")).reset_index()
    order = ["nonsense", "missense", "splice_region", "synonymous",
             "intronic", "5UTR", "3UTR", "inframe/indel", "start/stop_lost"]
    pool = pool.set_index("route_class").reindex(order).dropna(how="all").reset_index()
    y = np.arange(len(pool))
    ax[2].barh(y - 0.2, pool.P2, 0.4, label="P/LP", color="#c0392b")
    ax[2].barh(y + 0.2, pool.B2, 0.4, label="B/LB", color="#2e86c1")
    ax[2].set_yticks(y); ax[2].set_yticklabels(pool.route_class); ax[2].invert_yaxis()
    ax[2].set_xlabel("controls (n, pooled 2★+)")
    ax[2].set_title("c  Per-stratum controls"); ax[2].legend(frameon=False)

    fig.tight_layout()
    fig.savefig(paths.RESULTS / "control_summary.png", dpi=150, bbox_inches="tight")
    print("\nwrote control_summary.png")


def main() -> None:
    ctrl = load()
    join_diagnostics(ctrl)
    ct = concordance(ctrl)
    prelim_oddspath(ctrl)
    figure(ctrl, ct)


if __name__ == "__main__":
    main()
