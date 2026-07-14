#!/usr/bin/env python3
"""
figures.py — the one figure that carries the reliability-map finding.

Reads results/reliability_discrimination.csv and results/missed_demotions.csv
(produced by reliability.py) and draws a two-panel figure:

  A. Reliability vs distance to junction: AG splice AUROC by distance bin with
     bootstrap CI. The reliable band (CI-lower > 0.70) is shaded. Shows AG is
     reliable only at the junction edge and drops to chance mid-exon.
  B. Missed demotions by distance: fraction of measured demotions AG misses at
     90% specificity. Near-total mid-exon.

Output (results/): reliability_map.png

Run:  python figures.py
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analysis_common as C

RELIABLE = "#2166AC"
UNREL = "#B2182B"
GREY = "0.6"
ORDER = ["exonic_splice_region(≤3)", "4-10", "11-25", "26-50", "51-100", ">100"]
LABELS = ["≤3\n(edge)", "4-10", "11-25", "26-50", "51-100", ">100\n(mid-exon)"]


def main():
    disc = pd.read_csv(os.path.join(C.RESULTS_DIR, "reliability_discrimination.csv"))
    miss = pd.read_csv(os.path.join(C.RESULTS_DIR, "missed_demotions.csv"))
    d = disc[disc.ag_score == "ag_splice_mag"]

    # assemble the distance series: edge stratum + fine exon-core bins
    def row(stype, sval):
        r = d[(d.stratum_type == stype) & (d.stratum_value == sval)]
        return r.iloc[0] if len(r) else None
    series = [("dist_bin", "exonic_splice_region(≤3)")] + [("fine_bin", b) for b in ORDER[1:]]
    auroc = [row(t, v).auroc for t, v in series]
    lo = [row(t, v).auroc_lo for t, v in series]
    hi = [row(t, v).auroc_hi for t, v in series]
    colors = [RELIABLE if l > C.RELIABLE_CI_LO else UNREL for l in lo]

    m = miss[miss.ag_score.isna() if "ag_score" in miss else np.ones(len(miss), bool)] if False else miss
    def missrate(stype, sval):
        r = miss[(miss.stratum_type == stype) & (miss.stratum_value == sval)]
        return r.iloc[0].missed_rate_90spec if len(r) else np.nan
    mr = [missrate("dist_bin", "exonic_splice_region(≤3)")] + [missrate("fine_bin", b) for b in ORDER[1:]]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(10.5, 4.2))
    x = np.arange(len(series))

    # Panel A
    axA.axhspan(C.RELIABLE_CI_LO, 1.0, color=RELIABLE, alpha=0.06, zorder=0)
    axA.errorbar(x, auroc, yerr=[np.array(auroc) - np.array(lo), np.array(hi) - np.array(auroc)],
                 fmt="none", ecolor=GREY, capsize=3, lw=1, zorder=2)
    axA.scatter(x, auroc, c=colors, s=48, zorder=3, edgecolor="w", linewidth=0.6)
    axA.axhline(0.5, color="k", lw=0.8, ls=":")
    axA.axhline(C.RELIABLE_CI_LO, color=RELIABLE, lw=0.9, ls="--")
    axA.set_xticks(x); axA.set_xticklabels(LABELS, fontsize=7.5)
    axA.set_xlabel("distance to junction (nt)")
    axA.set_ylabel("AG splice AUROC vs measured mRNA-drop")
    axA.set_ylim(0.45, 0.9)
    axA.set_title("A  Reliability collapses with distance", loc="left", fontsize=10)
    axA.text(0.02, C.RELIABLE_CI_LO + 0.005, "reliable (CI-lo > 0.70)", color=RELIABLE, fontsize=6.8,
             transform=axA.get_yaxis_transform())

    # Panel B
    axB.bar(x, np.array(mr) * 100, color=colors, zorder=3, width=0.62)
    axB.set_xticks(x); axB.set_xticklabels(LABELS, fontsize=7.5)
    axB.set_xlabel("distance to junction (nt)")
    axB.set_ylabel("% of measured demotions AG misses (90% spec)")
    axB.set_ylim(0, 105)
    axB.set_title("B  Mid-exon demotions are invisible to AG", loc="left", fontsize=10)
    for xi, v in zip(x, mr):
        axB.annotate(f"{v*100:.0f}%", xy=(xi, v * 100), xytext=(0, 3),
                     textcoords="offset points", ha="center", fontsize=7)

    fig.tight_layout()
    out = os.path.join(C.RESULTS_DIR, "reliability_map.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
