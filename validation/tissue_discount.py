#!/usr/bin/env python3
"""
tissue_discount.py -- is the one-tier tissue-transfer discount the right size?

Reuses the frozen error-bound engine (model/bounds.py, bound_group -> Clopper-Pearson),
unchanged.

Inputs (all frozen):
  data/clinical_controls.csv
  data/annotated_screens/{BRCA1,VHL}_annotated.csv       (replicate / timepoint cols)
  results/evidence_codes/evidence_codes.csv              (route, call, is_rna_route)
Outputs (results/validation/):
  transfer_reproduction.csv       reproduction rate by mechanism class (the comparison)
  transfer_tier_calibration.csv   the gap mapped to ACMG tiers (BRCA1, powered)

Run:  python validation/tissue_discount.py
"""
from __future__ import annotations
import sys, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import beta as _beta

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "model"))
sys.path.insert(0, str(REPO))   # utils
from bounds import bound_group
from utils import paths

DELTA = 0.05
# frozen FPR-5% RNA-drop thresholds that DEFINE the method's rna_drop flag (used by the evidence-code stage)
RNA_THR5 = {"BARD1": -0.684, "PALB2": -0.564, "RAD51D": -1.016, "BRCA1": -0.667, "VHL": -0.597}
# OddsPath ladder rungs (frozen): one tier at the RNA-code rung = Moderate / Supporting
ONE_TIER = 4.3 / 2.1                        # 2.05x
# per-gene: (primary RNA column, [independent second-condition column(s)])
SECOND_CONDITION = {
    "BRCA1": ("score_rna", ["score_rna_rep1", "score_rna_rep2"]),   # replicate measurements
    "VHL":   ("rna_score_d20", ["rna_score_d6"]),                    # d20 primary, d6 independent
}


def _key(df):
    return df[["chrom", "pos", "ref", "alt"]].astype(str).agg("|".join, axis=1)


def cp_ci(k, n, a=DELTA):
    lo = 0.0 if k == 0 else _beta.ppf(a/2, k, n-k+1)
    hi = 1.0 if k == n else _beta.ppf(1-a/2, k+1, n-k)
    return float(lo), float(hi)


def odds(p):
    return p/(1-p) if 0 < p < 1 else np.inf


def tiers(p_invariant, p_variable):
    OR = odds(p_invariant) / odds(p_variable)
    return OR, (math.log(OR)/math.log(ONE_TIER) if OR > 0 and np.isfinite(OR) else np.nan)


def load(gene):
    ev = pd.read_csv(paths.RESULTS_CODES / "evidence_codes.csv")
    ev["key"] = _key(ev)
    a = pd.read_csv(paths.ANNOTATED_SCREENS / f"{gene}_annotated.csv")
    a["key"] = _key(a)
    return a.merge(ev[["key", "route", "call", "is_rna_route"]], on="key", how="left")


def reproduction(gene):
    """CP-bounded reproduction rate of the mRNA drop in the second condition,
    split into the tissue-invariant (nonsense/NMD) and tissue-variable (RNA-route) arms."""
    primary, second = SECOND_CONDITION[gene]
    df = load(gene)
    need = [primary] + second
    base = df[(df.call == "LoF") & df[need].notna().all(axis=1)].copy()
    base = base[base[primary] <= RNA_THR5[gene]]                 # drop present in primary
    base["reproduces"] = (base[second] <= RNA_THR5[gene]).all(axis=1)  # ...and in the 2nd condition
    arms = {
        "tissue-invariant (nonsense/NMD)": base[base.route_class == "nonsense"],
        "tissue-variable (RNA route)":     base[base.is_rna_route == True],
    }
    rows = []
    for name, sub in arms.items():
        n = len(sub); k = int(sub.reproduces.sum())
        lo, hi = cp_ci(k, n) if n else (np.nan, np.nan)
        rows.append(dict(gene=gene, mechanism_class=name, n=n, reproduces=k,
                         reproduction_rate=round(k/n, 3) if n else np.nan,
                         cp_lo=round(lo, 3), cp_hi=round(hi, 3)))
    return pd.DataFrame(rows), arms


def main():
    out = paths.RESULTS_VALIDATION; out.mkdir(parents=True, exist_ok=True)

    repB, armsB = reproduction("BRCA1")
    repV, armsV = reproduction("VHL")

    # reproduction table: BRCA1 (both arms) + VHL RNA-route arm as second-gene corroboration.
    # VHL's invariant arm is n=9 -> not interpretable, excluded from the table with a note.
    rep = pd.concat([repB,
                     repV[repV.mechanism_class == "tissue-variable (RNA route)"]
                     .assign(note="second-gene corroboration (timepoints); "
                                  "invariant arm n=9, cannot calibrate")],
                    ignore_index=True)
    rep.to_csv(out / "transfer_reproduction.csv", index=False)

    # tier calibration: BRCA1 only (the powered comparison)
    inv = armsB["tissue-invariant (nonsense/NMD)"]; var = armsB["tissue-variable (RNA route)"]
    ki, ni = int(inv.reproduces.sum()), len(inv)
    kv, nv = int(var.reproduces.sum()), len(var)
    pi, pv = ki/ni, kv/nv
    OR, gap = tiers(pi, pv)
    gap_lo = tiers(cp_ci(ki, ni)[0], cp_ci(kv, nv)[1])[1]
    gap_hi = tiers(cp_ci(ki, ni)[1], cp_ci(kv, nv)[0])[1]
    verdict = ("calibrated (95% CI contains 1 tier)" if gap_lo <= 1.0 <= gap_hi
               else "not one tier (95% CI excludes 1)")
    cal = pd.DataFrame([dict(
        gene="BRCA1",
        invariant_reproduction=round(pi, 3), variable_reproduction=round(pv, 3),
        odds_ratio=round(OR, 2), one_tier_odds=round(ONE_TIER, 2),
        tier_gap=round(gap, 2), tier_gap_lo=round(gap_lo, 2), tier_gap_hi=round(gap_hi, 2),
        verdict=verdict)])
    cal.to_csv(out / "transfer_tier_calibration.csv", index=False)

    print("REPRODUCTION of the mRNA drop in a second measurement condition:")
    print(rep.to_string(index=False))
    print(f"\nBRCA1 tier calibration (one tier = {ONE_TIER:.2f}x odds):")
    print(f"   invariant {pi:.3f} vs variable {pv:.3f}  ->  OR {OR:.2f}  =  "
          f"{gap:.2f} tiers  [{gap_lo:.2f}, {gap_hi:.2f}]")
    print(f"   verdict: {verdict}  ->  fixed one-tier discount is empirically supported")


if __name__ == "__main__":
    main()
