"""
build_strength.py -- driver: clinical OddsPath strength with a finite-sample error bound.

Reads the clinical controls and produces the error-bounded ACMG/AMP functional
strengths, using the Clopper-Pearson error-bound engine (bounds.py) unchanged and
the OddsPath map (oddspath.py). All deliverables are written to results/bounds/.


Outputs:
  oddspath_calls_per_class.csv   -- PRIMARY: gene x class x direction error-bounded
                                    OddsPath with lower bound, tier, abstention,
                                    and apply-variant coverage.
  oddspath_per_gene.csv          -- whole-gene summary (mixes classes; supersedes
                                    the preliminary (unbounded) per-gene table).
  oddspath_robustness.csv        -- full vs fit vs validate for non-abstaining
                                    per-class groups.
  fig_step2_oddspath.png         -- error-bounded OddsPath vs the ACMG ladder.
  STEP2_REPORT.md is written by hand from these tables.
"""

from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# make the module resolvable whether launched as `python step2_oddspath.py`
sys.path.insert(0, str(Path(__file__).resolve().parent))          # model/ (strength.py)
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))       # repo root (utils)
from strength import (oddspath_grid, oddspath_group, tier_from_oddspath,
                      ODDS_SUPPORTING, ODDS_MODERATE, ODDS_STRONG)

from utils import paths
HERE = paths.RESULTS_STRENGTH
CONTROLS = paths.CLINICAL_CONTROLS


def load_calibration() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (calibration controls, all assayed) with chrom as str."""
    ctrl = pd.read_csv(CONTROLS)
    ctrl["chrom"] = ctrl.chrom.astype(str)
    cal = ctrl[ctrl.role == "calibration"].copy()
    return cal, ctrl


def per_class_table(cal: pd.DataFrame, ctrl: pd.DataFrame) -> pd.DataFrame:
    grid = oddspath_grid(cal)

    # apply-variant coverage: how many reclassification targets sit in each
    # gene x class and would receive the code when called in that direction.
    apply = ctrl[ctrl.role == "apply"]
    cov = (apply.groupby(["gene", "route_class"])
           .agg(apply_total=("call", "size"),
                apply_LoF=("call", lambda s: int((s == "LoF").sum())),
                apply_Normal=("call", lambda s: int((s == "Normal").sum())))
           .reset_index())
    grid = grid.merge(cov, on=["gene", "route_class"], how="left")
    # variants the code would actually stamp: LoF for PS3, Normal for BS3
    grid["apply_receiving"] = np.where(
        grid.direction == "PS3", grid.apply_LoF, grid.apply_Normal)
    return grid


def per_gene_table(cal: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for g, sub in cal.groupby("gene"):
        for d in ("PS3", "BS3"):
            rec = dict(gene=g)
            rec.update(oddspath_group(sub, d))
            rows.append(rec)
    lead = ["gene", "direction", "n_P", "n_B", "oddspath_point",
            "oddspath_lo", "tier", "abstain", "abstain_reason"]
    df = pd.DataFrame(rows)
    return df[lead + [c for c in df.columns if c not in lead]]


def robustness_table(cal: pd.DataFrame, per_class: pd.DataFrame) -> pd.DataFrame:
    keep = per_class[~per_class.abstain][["gene", "route_class", "direction"]]
    rows = []
    for _, r in keep.iterrows():
        g, cl, d = r.gene, r.route_class, r.direction
        sub = cal[(cal.gene == g) & (cal.route_class == cl)]
        rec = dict(gene=g, route_class=cl, direction=d)
        for half in ("full", "fit", "validate"):
            s = sub if half == "full" else sub[sub.calib_fold == half]
            res = oddspath_group(s, d)
            rec[f"{half}_oddspath_lo"] = res["oddspath_lo"]
            rec[f"{half}_tier"] = res["tier"]
            rec[f"{half}_nP"] = res["n_P"]
            rec[f"{half}_nB"] = res["n_B"]
        rows.append(rec)
    return pd.DataFrame(rows)


GENE_ORDER = ["BARD1", "BRCA1", "PALB2", "RAD51D", "VHL"]
GCOL = dict(zip(GENE_ORDER,
                ["#c0392b", "#2e86c1", "#27ae60", "#8e44ad", "#e67e22"]))
XMAX = 700.0


def _forest(ax, rows, ylabels):
    """One forest panel: rows = [((gene, key), [(gene, dir, point, lo), ...]), ...].
    Filled marker = error-bounded lower bound (reported); open marker + connector =
    point estimate; arrow to the panel edge when the point estimate is infinite.
    Marker shape: circle = PS3, square = BS3. Bands = ACMG strength ladder."""
    bands = [(1, ODDS_SUPPORTING, "#f2f2f2"),
             (ODDS_SUPPORTING, ODDS_MODERATE, "#fde5cf"),
             (ODDS_MODERATE, ODDS_STRONG, "#f9c784"),
             (ODDS_STRONG, XMAX, "#f0995e")]
    for lo, hi, c in bands:
        ax.axvspan(lo, hi, color=c, zorder=0)
    for thr in (ODDS_SUPPORTING, ODDS_MODERATE, ODDS_STRONG):
        ax.axvline(thr, color="k", lw=0.6, ls=":", zorder=1)
    for y, (key, recs) in enumerate(rows):
        for g, d, pt, lo in recs:
            mk = "o" if d == "PS3" else "s"
            ax.scatter([lo], [y], marker=mk, s=70, color=GCOL[g],
                       edgecolor="k", linewidth=0.6, zorder=4)
            if np.isfinite(pt) and pt < XMAX:
                ax.plot([lo, pt], [y, y], color=GCOL[g], lw=1.1,
                        alpha=0.55, zorder=3)
                ax.scatter([pt], [y], marker=mk, s=34, facecolor="white",
                           edgecolor=GCOL[g], linewidth=1.0, zorder=3)
            else:
                ax.annotate("", xy=(XMAX * 0.96, y), xytext=(lo, y),
                            arrowprops=dict(arrowstyle="-|>", color=GCOL[g],
                                            lw=1.1, alpha=0.55), zorder=3)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(ylabels)
    for t, (key, _) in zip(ax.get_yticklabels(), rows):
        t.set_color(GCOL[key[0]])
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(1, XMAX)
    ax.set_xticks([1, ODDS_SUPPORTING, ODDS_MODERATE, ODDS_STRONG, 100, 700])
    ax.set_xticklabels(["1", "2.1", "4.3", "18.7", "100", "700"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def figure(per_class: pd.DataFrame, per_gene: pd.DataFrame) -> None:
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 9,
                         "axes.labelsize": 9, "xtick.labelsize": 7,
                         "ytick.labelsize": 8, "axes.titlelocation": "left"})

    # panel a: per-class error-bounded codes, one row per gene x class, both dirs
    codes = per_class[~per_class.abstain]
    gxc = list(dict.fromkeys(zip(codes.gene, codes.route_class)))
    rows_a, lab_a = [], []
    for g, cl in gxc:
        recs = []
        for d in ("PS3", "BS3"):
            r = codes[(codes.gene == g) & (codes.route_class == cl) &
                      (codes.direction == d)]
            if len(r):
                recs.append((g, d, r.oddspath_point.iloc[0], r.oddspath_lo.iloc[0]))
        rows_a.append(((g, cl), recs))
        lab_a.append(f"{g}  {cl.replace('_', '-')}")

    # panel b: whole-gene summary
    rows_b, lab_b = [], []
    for g in GENE_ORDER:
        recs = [(g, d,
                 per_gene[(per_gene.gene == g) & (per_gene.direction == d)].oddspath_point.iloc[0],
                 per_gene[(per_gene.gene == g) & (per_gene.direction == d)].oddspath_lo.iloc[0])
                for d in ("PS3", "BS3")]
        rows_b.append(((g, "whole"), recs))
        lab_b.append(g)

    fig = plt.figure(figsize=(9.2, 7.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[len(rows_a), len(rows_b)],
                          hspace=0.30)
    axa, axb = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    _forest(axa, rows_a, lab_a)
    _forest(axb, rows_b, lab_b)
    axa.set_title("a   Per-class error-bounded codes — protein-blind splice-region "
                  "earns a code in every gene (Moderate in 4, Supporting in VHL), no pooling")
    axb.set_title("b   Whole-gene summary (mixes classes) — supersedes the "
                  "preliminary (unbounded) estimates")
    axb.set_xlabel("error-bounded OddsPath — conservative lower bound   "
                   "(larger = stronger evidence)")
    for x, name in [(np.sqrt(1 * ODDS_SUPPORTING), "below"),
                    (np.sqrt(ODDS_SUPPORTING * ODDS_MODERATE), "Supporting"),
                    (np.sqrt(ODDS_MODERATE * ODDS_STRONG), "Moderate"),
                    (np.sqrt(ODDS_STRONG * XMAX), "Strong")]:
        axb.annotate(name, xy=(x, 0), xytext=(x, -1.35),
                     textcoords=("data", "offset fontsize"), ha="center",
                     va="top", fontsize=7, color="#666", annotation_clip=False)

    from matplotlib.lines import Line2D
    leg = [Line2D([], [], marker="o", color="none", markerfacecolor="#555",
                  markeredgecolor="k", markersize=8, label="PS3 (pathogenic)"),
           Line2D([], [], marker="s", color="none", markerfacecolor="#555",
                  markeredgecolor="k", markersize=8, label="BS3 (benign)"),
           Line2D([], [], marker="o", color="none", markerfacecolor="white",
                  markeredgecolor="#555", markersize=7, label="point estimate"),
           Line2D([], [], marker="o", color="none", markerfacecolor="#555",
                  markeredgecolor="k", markersize=8, label="error-bounded bound")]
    axa.legend(handles=leg, loc="lower right", frameon=False, fontsize=7,
               handletextpad=0.4)

    fig.savefig(HERE / "fig_step2_oddspath.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    cal, ctrl = load_calibration()
    per_class = per_class_table(cal, ctrl)
    per_gene = per_gene_table(cal)
    robust = robustness_table(cal, per_class)

    per_class.to_csv(HERE / "oddspath_calls_per_class.csv", index=False)
    per_gene.to_csv(HERE / "oddspath_per_gene.csv", index=False)
    robust.to_csv(HERE / "oddspath_robustness.csv", index=False)
    figure(per_class, per_gene)

    n_codes = int((~per_class.abstain).sum())
    print(f"per-class groups: {len(per_class)}  error-bounded (non-abstain): {n_codes}")
    print(f"whole-gene rows : {len(per_gene)}")
    print(f"robustness rows : {len(robust)}")
    print("wrote oddspath_calls_per_class.csv, oddspath_per_gene.csv, "
          "oddspath_robustness.csv, fig_step2_oddspath.png")


if __name__ == "__main__":
    main()
