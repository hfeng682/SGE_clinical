#!/usr/bin/env python3
"""
bound_holds.py -- does the finite-sample error bound actually hold on held-out controls?

Two held-out tests, both reusing the FROZEN error-bound engine (model/bounds.py)
unchanged:

  1. WITHIN-GENE (how the method is used). Per gene x stratum, set the call-error
     ceiling (Clopper-Pearson, position-clustered) on the calibration FIT controls,
     then measure the observed error on the held-out calibration VALIDATE controls.
     covered = observed <= ceiling. Every gene is bounded on its own controls
     (the no-pooling design). Pre-registered target: empirical coverage >= 1-delta.

  2. LEAVE-ONE-GENE-OUT (negative control). For each held-out gene, set the ceiling
     from the OTHER four genes' pooled controls per stratum, then test whether it
     covers the held-out gene. Because it deliberately pools across genes, it is a
     direct test of the no-pooling design: if pooling covered, genes could borrow
     strength; it does not, which is why the method never pools.

The two error rates bounded per stratum:
  false-loss   = benign 2* control called LoF        -> controls the PS3 claim
  false-normal = pathogenic 2* control called Normal -> controls the BS3 claim

Inputs (all frozen, nothing refit):
  data/clinical_controls.csv        role, calib_fold, call, clin_label, pos, ...
Outputs (results/validation/):
  within_gene_coverage.csv, logo_cross_gene_coverage.csv, logo_accumulation_curve.csv,
  within_gene_heldout_agreement.csv, within_gene_vus_graded.csv

Run:  python validation/bound_holds.py
"""
from __future__ import annotations
import sys
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "model"))
sys.path.insert(0, str(REPO))   # utils
from bounds import bound_group
from utils import paths

DELTA = 0.05
GENES = ["BARD1", "BRCA1", "PALB2", "RAD51D", "VHL"]
MIN_FIT, MIN_VAL, MIN_DONOR, MIN_TEST = 10, 5, 10, 5   # power floors


def _err_indicator(sub: pd.DataFrame, side: str) -> np.ndarray:
    """1 where the control is mis-called: benign->LoF (FPR) / pathogenic->Normal (FNR)."""
    err_call = "LoF" if side == "FPR" else "Normal"
    return (sub.call == err_call).astype(float).values


def _controls(df, gene, stratum, side):
    lab = "B/LB" if side == "FPR" else "P/LP"
    return df[(df.gene == gene) & (df.route_class == stratum) & (df.clin_label == lab)]


def within_gene(cal: pd.DataFrame) -> pd.DataFrame:
    fit, val = cal[cal.calib_fold == "fit"], cal[cal.calib_fold == "validate"]
    strata = sorted(cal.route_class.dropna().unique())
    rows = []
    for g in GENES:
        for st in strata:
            for side in ("FPR", "FNR"):
                f, v = _controls(fit, g, st, side), _controls(val, g, st, side)
                if len(f) < MIN_FIT or len(v) < MIN_VAL:
                    continue
                cg = bound_group(_err_indicator(f, side), f.pos.values, DELTA)
                obs = float(_err_indicator(v, side).mean())
                rows.append(dict(gene=g, stratum=st, side=side,
                                 direction={"FPR": "PS3", "FNR": "BS3"}[side],
                                 n_fit=len(f), n_val=len(v),
                                 fit_ceiling=round(cg["cp_ceiling"], 4),
                                 val_observed=round(obs, 4),
                                 covered=bool(obs <= cg["cp_ceiling"])))
    return pd.DataFrame(rows)


def leave_one_gene_out(cal: pd.DataFrame) -> pd.DataFrame:
    strata = sorted(cal.route_class.dropna().unique())
    rows = []
    for held in GENES:
        donors, test = cal[cal.gene != held], cal[cal.gene == held]
        for st in strata:
            for side in ("FPR", "FNR"):
                lab = "B/LB" if side == "FPR" else "P/LP"
                d = donors[(donors.route_class == st) & (donors.clin_label == lab)]
                t = test[(test.route_class == st) & (test.clin_label == lab)]
                if len(d) < MIN_DONOR or len(t) < MIN_TEST:
                    continue
                cg = bound_group(_err_indicator(d, side), d.pos.values, DELTA)
                if cg["abstain"]:
                    continue
                obs = float(_err_indicator(t, side).mean())
                rows.append(dict(held_out_gene=held, stratum=st,
                                 direction={"FPR": "PS3", "FNR": "BS3"}[side],
                                 n_donor=len(d), n_test=len(t),
                                 donor_ceiling=round(cg["cp_ceiling"], 4),
                                 heldout_observed=round(obs, 4),
                                 covered=bool(obs <= cg["cp_ceiling"])))
    return pd.DataFrame(rows)


def heldout_agreement(cal: pd.DataFrame) -> pd.DataFrame:
    """Scientific payload A: do the held-out (validate-fold) calls agree with ClinVar truth?
    These controls never informed any threshold, so this is honest clinical accuracy."""
    val = cal[cal.calib_fold == "validate"]
    rows = []
    for g in GENES + ["(pooled)"]:
        d = val if g == "(pooled)" else val[val.gene == g]
        B, P = d[d.clin_label == "B/LB"], d[d.clin_label == "P/LP"]
        rows.append(dict(gene=g, n_benign=len(B),
                         benign_to_Normal=round((B.call == "Normal").mean(), 3) if len(B) else np.nan,
                         n_path=len(P),
                         path_to_LoF=round((P.call == "LoF").mean(), 3) if len(P) else np.nan))
    return pd.DataFrame(rows)


def vus_graded(cal: pd.DataFrame, controls: pd.DataFrame) -> pd.DataFrame:
    """Scientific payload B: how many real held-out (apply-role) variants receive a error-bounded
    ACMG code from the 11 error-bounded strata? Reads the frozen error-bounded-code table."""
    codes = pd.read_csv(paths.RESULTS_STRENGTH / "oddspath_calls_per_class.csv")
    cert = codes[codes.tier.isin(["Supporting", "Moderate", "Strong"])]
    apply = controls[controls.role == "apply"]
    rows = []
    for (g, st), sub in cert.groupby(["gene", "route_class"]):
        dirs = sub.direction.tolist()
        ap = apply[(apply.gene == g) & (apply.route_class == st)]
        n_lof, n_norm = int((ap.call == "LoF").sum()), int((ap.call == "Normal").sum())
        rows.append(dict(gene=g, stratum=st, code_earned="+".join(sorted(set(dirs))),
                         apply_total=len(ap),
                         graded_VUS=(n_lof if "PS3" in dirs else 0) + (n_norm if "BS3" in dirs else 0),
                         as_LoF_PS3=n_lof if "PS3" in dirs else 0,
                         as_Normal_BS3=n_norm if "BS3" in dirs else 0))
    return pd.DataFrame(rows)


def accumulation(cal: pd.DataFrame) -> pd.DataFrame:
    """Cross-gene coverage as a function of #donor genes (avg over donor subsets)."""
    strata = sorted(cal.route_class.dropna().unique())
    acc = {k: [] for k in range(1, len(GENES))}
    for held in GENES:
        others = [x for x in GENES if x != held]
        test = cal[cal.gene == held]
        for st in strata:
            for side in ("FPR", "FNR"):
                lab = "B/LB" if side == "FPR" else "P/LP"
                t = test[(test.route_class == st) & (test.clin_label == lab)]
                if len(t) < MIN_TEST:
                    continue
                obs = float(_err_indicator(t, side).mean())
                for k in range(1, len(GENES)):
                    for donor_set in combinations(others, k):
                        d = cal[(cal.gene.isin(donor_set)) &
                                (cal.route_class == st) & (cal.clin_label == lab)]
                        if len(d) < MIN_DONOR:
                            continue
                        cg = bound_group(_err_indicator(d, side), d.pos.values, DELTA)
                        if cg["abstain"]:
                            continue
                        acc[k].append(obs <= cg["cp_ceiling"])
    return pd.DataFrame([dict(n_donor_genes=k, coverage=round(np.mean(v), 4), n_evals=len(v))
                         for k, v in acc.items()])


def main() -> None:
    controls = pd.read_csv(paths.CLINICAL_CONTROLS)
    controls["chrom"] = controls.chrom.astype(str)
    cal = controls[controls.role == "calibration"].copy()
    out = paths.RESULTS_VALIDATION; out.mkdir(parents=True, exist_ok=True)

    wg = within_gene(cal)
    wg.to_csv(out / "within_gene_coverage.csv", index=False)
    logo = leave_one_gene_out(cal)
    logo.to_csv(out / "logo_cross_gene_coverage.csv", index=False)
    acc = accumulation(cal)
    acc.to_csv(out / "logo_accumulation_curve.csv", index=False)
    ag = heldout_agreement(cal)
    ag.to_csv(out / "within_gene_heldout_agreement.csv", index=False)
    vg = vus_graded(cal, controls)
    vg.to_csv(out / "within_gene_vus_graded.csv", index=False)

    print(f"WITHIN-GENE coverage : {wg.covered.mean():.3f} "
          f"({int(wg.covered.sum())}/{len(wg)})  target >= {1-DELTA}  "
          f"[{'PASS' if wg.covered.mean() >= 1-DELTA else 'FAIL'}]")
    print(f"LOGO cross-gene      : {logo.covered.mean():.3f} "
          f"({int(logo.covered.sum())}/{len(logo)})  target >= {1-DELTA}  "
          f"[{'PASS' if logo.covered.mean() >= 1-DELTA else 'FAIL'}]")
    print("accumulation (coverage vs #donor genes):")
    for _, r in acc.iterrows():
        print(f"   {int(r.n_donor_genes)} gene(s): {r.coverage:.3f}  (n={int(r.n_evals)})")
    pooled = ag[ag.gene == "(pooled)"].iloc[0]
    print(f"held-out agreement   : path->LoF {pooled.path_to_LoF}, "
          f"benign->Normal {pooled.benign_to_Normal} (validate fold, vs ClinVar)")
    print(f"VUS graded by 11 codes: {int(vg.graded_VUS.sum())} real held-out variants "
          f"({int(vg[vg.stratum=='splice_region'].graded_VUS.sum())} protein-blind splice)")


if __name__ == "__main__":
    main()
