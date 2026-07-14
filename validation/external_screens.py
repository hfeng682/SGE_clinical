#!/usr/bin/env python3
"""
external_screens.py -- does the survival axis transfer to screens the method never saw?

The method has two axes:
  Axis 1  survival-score three-class call -> Clopper-Pearson ceiling -> OddsPath
          -> ACMG strength, with the distribution-free finite-sample bound.
  Axis 2  mechanism route (protein / RNA) + tissue-transfer discount, which REQUIRES
          a second (mRNA-abundance) score.

No genuinely external *two-score* SGE screen exists yet (every mRNA-readout SGE is
HAP1-lineage; see the two-score landscape table (data/sge_two_score_landscape.csv)). But Axis 1 needs only a
survival/function score, so it CAN be exercised on external one-score screens. This
script does exactly that on two screens the method never saw, taking each lab's OWN
published three-class call as input (the method consumes an investigator call; it does
not re-derive one):

  RAD51C  urn:mavedb:00000673-0-1  (Olvera-Leon et al., Cell 2024)
          fast/slow depleted -> LoF; unchanged/enriched -> Normal
  DDX3X   urn:mavedb:00000658-0-1  (X-linked neurodevelopmental disorder)
          normal -> Normal; abnormal+High -> LoF; abnormal+Intermediate -> Uncertain

Only Axis 1 is reported. The mechanism/tissue axis is deliberately absent -- these
screens have no second score -- and that limitation is the point: it is what motivates
two-score screening.

Inputs staged in validation/external_data/ (see external_data/PROVENANCE.md):
  external_data/annotated/<GENE>_annotated.csv   VEP-annotated screen (route_class, pos, ...)
  external_data/raw/<GENE>_scores.csv            MaveDB score set (investigator call columns)
  external_data/controls/external_clinical_controls.csv   pre-built controls (this script rebuilds it)
Outputs (results/validation/):
  external_axis1_oddspath.csv, external_clinical_agreement.csv,
  external_per_stratum_concordance.csv, external_vus_graded.csv

Run:  python validation/external_screens.py
"""
from __future__ import annotations
import sys, gzip, json
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "model"))
sys.path.insert(0, str(REPO))   # utils
from strength import oddspath_grid
from utils import paths

VCF = paths.CLINVAR_VCF
PAD = 200
EXTERNAL = {  # gene -> (chrom, MaveDB urn)
    "RAD51C": ("17", "urn:mavedb:00000673-0-1"),
    "DDX3X":  ("X",  "urn:mavedb:00000658-0-1"),
}

# ClinVar review status -> gold stars, and CLNSIG -> truth bucket (identical to the clinical-controls stage)
STARS = {"practice_guideline": 4, "reviewed_by_expert_panel": 3,
         "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "criteria_provided,_single_submitter": 1,
         "criteria_provided,_conflicting_classifications": 1,
         "criteria_provided,_conflicting_interpretations": 1}
P_SET = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
B_SET = {"Benign", "Likely_benign", "Benign/Likely_benign"}
LAB_RANK = {"P/LP": 0, "B/LB": 1, "conflicting": 2, "VUS": 3, "other": 4}


def clnsig_bucket(s):
    if s in P_SET: return "P/LP"
    if s in B_SET: return "B/LB"
    if s == "Uncertain_significance": return "VUS"
    if isinstance(s, str) and s.startswith("Conflicting"): return "conflicting"
    return "other"


def investigator_call(gene: str, raw: pd.DataFrame) -> pd.Series:
    """Map each lab's PUBLISHED functional class onto {Normal, Uncertain, LoF}."""
    if gene == "RAD51C":
        m = {"fast depleted": "LoF", "slow depleted": "LoF",
             "unchanged": "Normal", "enriched": "Normal"}
        return raw["functional_classification"].map(m)
    if gene == "DDX3X":
        def call(pred, conf):
            if pred == "normal": return "Normal"
            conf = None if (pd.isna(conf) or conf == "NA") else conf
            return {"High_confidence": "LoF",
                    "Intermediate_confidence": "Uncertain"}.get(conf, "Uncertain")
        return pd.Series([call(p, c) for p, c in zip(
            raw["SGE_prediction_of_variant_function_in_NDD_context"],
            raw["Confidence_of_functionally_abnormal_variant_prediction"])], index=raw.index)
    raise ValueError(gene)


def scan_clinvar(regions: dict) -> pd.DataFrame:
    windows = {}
    for g, (c, lo, hi) in regions.items():
        windows.setdefault(c, []).append((g, lo - PAD, hi + PAD))
    chroms = set(windows)
    recs = []
    with gzip.open(VCF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom = line[:line.find("\t")]
            if chrom not in chroms:
                continue
            p = line.rstrip("\n").split("\t"); pos = int(p[1])
            gene = next((gn for gn, lo, hi in windows[chrom] if lo <= pos <= hi), None)
            if gene is None:
                continue
            info = {k: v for kv in p[7].split(";") if "=" in kv for k, v in [kv.split("=", 1)]}
            recs.append(dict(gene=gene, chrom=chrom, pos=pos, ref=p[3], alt=p[4],
                             clinvar_id=p[2], CLNSIG=info.get("CLNSIG", ""),
                             gold_stars=STARS.get(info.get("CLNREVSTAT", ""), 0)))
    cv = pd.DataFrame(recs)
    cv["clin_label"] = cv.CLNSIG.map(clnsig_bucket)
    return cv


def build_controls(gene, screen, cv, gnomad):
    a = screen.copy(); a["chrom"] = a.chrom.astype(str); a["gene"] = gene
    c = cv[cv.gene == gene].assign(_lab=lambda d: d.clin_label.map(LAB_RANK))
    best = (c.sort_values(["gold_stars", "_lab"], ascending=[False, True])
              .drop_duplicates(["chrom", "pos", "ref", "alt"]))
    key = ["chrom", "pos", "ref", "alt"]
    j = a.merge(best[key + ["clin_label", "gold_stars", "clinvar_id"]], on=key, how="left")
    if gnomad is not None and len(gnomad):
        j = j.merge(gnomad, on=key, how="left")
    j["gold_stars"] = j.gold_stars.fillna(0).astype(int)
    j["role"] = np.where((j.gold_stars >= 2) & j.clin_label.isin(["P/LP", "B/LB"]),
                         "calibration", "apply")
    return j


def main() -> None:
    out = paths.RESULTS_VALIDATION; out.mkdir(parents=True, exist_ok=True)

    # 1. load annotated screens + attach investigator call
    screens, regions = {}, {}
    for g in EXTERNAL:
        ann = pd.read_csv(paths.REPO / "validation" / "external_data" / "annotated" / f"{g}_annotated.csv")
        ann["chrom"] = ann.chrom.astype(str)
        raw = pd.read_csv(paths.REPO / "validation" / "external_data" / "raw" / f"{g}_scores.csv")
        # annotated files already carry the investigator columns (passed through by VEP front-end)
        src = ann if any(col in ann.columns for col in
                         ("functional_classification",
                          "SGE_prediction_of_variant_function_in_NDD_context")) else \
              ann.merge(raw, on="hgvs_nt", how="left")
        ann["call"] = investigator_call(g, src)
        screens[g] = ann
        regions[g] = (str(ann.chrom.mode().iloc[0]), int(ann.pos.min()), int(ann.pos.max()))

    # 2. ClinVar controls from the local VCF (gnomAD AF optional -- pre-built controls carries it)
    cv = scan_clinvar(regions)

    # 3. per-gene controls + Axis 1 error bounding via the frozen engine
    certs, agree, conc, yields = [], [], [], []
    for g in EXTERNAL:
        ctrl = build_controls(g, screens[g], cv, gnomad=None)
        cal = ctrl[(ctrl.role == "calibration") & ctrl.call.notna()
                  & ctrl.clin_label.isin(["P/LP", "B/LB"])]
        grid = oddspath_grid(cal); grid["gene"] = g
        certs.append(grid)

        # clinical agreement vs independent 2* ClinVar
        B, P = cal[cal.clin_label == "B/LB"], cal[cal.clin_label == "P/LP"]
        agree.append(dict(gene=g, n_benign=len(B),
                          benign_to_Normal=int((B.call == "Normal").sum()),
                          benign_agree=round((B.call == "Normal").mean(), 3) if len(B) else np.nan,
                          n_path=len(P), path_to_LoF=int((P.call == "LoF").sum()),
                          path_agree=round((P.call == "LoF").mean(), 3) if len(P) else np.nan))
        for st, sub in cal.groupby("route_class"):
            Bs, Ps = sub[sub.clin_label == "B/LB"], sub[sub.clin_label == "P/LP"]
            conc.append(dict(gene=g, stratum=st, n_P=len(Ps), n_B=len(Bs),
                             path_to_LoF=f"{int((Ps.call=='LoF').sum())}/{len(Ps)}" if len(Ps) else "-",
                             benign_to_Normal=f"{int((Bs.call=='Normal').sum())}/{len(Bs)}" if len(Bs) else "-"))
        # VUS graded by error-bounded codes
        cert = grid[~grid.abstain]
        for st, sub in cert.groupby("route_class"):
            dirs = sub.direction.tolist()
            ap = ctrl[(ctrl.role == "apply") & (ctrl.route_class == st)]
            n_lof, n_norm = int((ap.call == "LoF").sum()), int((ap.call == "Normal").sum())
            yields.append(dict(gene=g, stratum=st, code_earned="+".join(sorted(dirs)),
                               apply_total=len(ap),
                               graded_VUS=(n_lof if "PS3" in dirs else 0) + (n_norm if "BS3" in dirs else 0)))

    ext_grid = pd.concat(certs, ignore_index=True)
    ext_grid.to_csv(out / "external_axis1_oddspath.csv", index=False)
    pd.DataFrame(agree).to_csv(out / "external_clinical_agreement.csv", index=False)
    pd.DataFrame(conc).to_csv(out / "external_per_stratum_concordance.csv", index=False)
    pd.DataFrame(yields).to_csv(out / "external_vus_graded.csv", index=False)

    for g in EXTERNAL:
        cert = ext_grid[(ext_grid.gene == g) & (~ext_grid.abstain)]
        a = next(r for r in agree if r["gene"] == g)
        print(f"{g}: {len(cert)} error-bounded code(s); "
              f"benign->Normal {a['benign_agree']}, path->LoF {a['path_agree']} "
              f"(vs independent 2* ClinVar)")


if __name__ == "__main__":
    main()
