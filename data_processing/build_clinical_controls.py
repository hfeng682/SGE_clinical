#!/usr/bin/env python3
"""

Build the clinical control set for the two-score SGE project.

This script performs the DATA PROCESSING for the clinical-controls stage. It assumes the raw sources
have already been downloaded/cached into this directory (see STEP0_REPORT.md,
"How the controls was built"):

    controls/clinvar.vcf.gz          ClinVar GRCh38 VCF (downloaded from NCBI)
    controls/gnomad_r4_<GENE>.json   gnomAD r4 region pulls (from the variants connector)

and that the annotated screens and the pooled three-class calls exist upstream:

    out/<GENE>_annotated.csv       VEP-annotated screens (from sge_pipeline/)
    analysis/pooled_labeled.csv    the shared {Normal/Uncertain/LoF} call per variant

Run from inside the controls/ directory:

    python build_clinical_controls.py

Outputs (written next to this script):
    gene_regions.csv
    clinvar_raw_regions.csv
    clinvar_labeled_regions.csv
    assay_clinvar_joined.csv
    clinical_controls.csv
    overlap_per_gene.csv
    overlap_per_stratum.csv
    split_manifest.json
"""
from __future__ import annotations
import gzip, json, hashlib, datetime
from pathlib import Path
import pandas as pd
import numpy as np

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils import paths                              # single source of truth for paths

DERIV  = paths.CLINVAR_GNOMAD                        # derivation intermediates + control-count tables
OUT    = paths.ANNOTATED_SCREENS                     # VEP-annotated screens
POOLED = paths.POOLED                                # shared three-class call per variant

GENES  = ["BARD1", "PALB2", "RAD51D", "BRCA1", "VHL"]
PAD    = 200            # bp padding added to each gene's assayed span for the VCF scan
SEED   = 20260710       # frozen split seed

# ClinVar review status -> gold stars (0-4).
STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    # everything else (no assertion criteria, blank, etc.) -> 0
}
# CLNSIG aggregate germline classification -> truth bucket.
P_SET = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
B_SET = {"Benign", "Likely_benign", "Benign/Likely_benign"}


def clnsig_bucket(sig: str) -> str:
    if sig in P_SET:                                     return "P/LP"
    if sig in B_SET:                                     return "B/LB"
    if sig == "Uncertain_significance":                  return "VUS"
    if sig.startswith("Conflicting"):                    return "conflicting"
    return "other"


def parse_info(info: str) -> dict:
    out = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


# 1. Gene coordinate windows (from the assayed variants themselves)
def gene_regions() -> pd.DataFrame:
    rows = []
    for g in GENES:
        d = pd.read_csv(OUT / f"{g}_annotated.csv", usecols=["chrom", "pos"])
        rows.append(dict(gene=g, chrom=str(d.chrom.iloc[0]),
                         n=len(d), pos_min=int(d.pos.min()), pos_max=int(d.pos.max())))
    reg = pd.DataFrame(rows)
    reg.to_csv(DERIV / "gene_regions.csv", index=False)
    return reg


# 2. Scan the ClinVar VCF once, keeping only records inside a gene window
def parse_clinvar(reg: pd.DataFrame) -> pd.DataFrame:
    # chrom -> list of (gene, lo, hi)
    windows: dict[str, list] = {}
    for _, r in reg.iterrows():
        windows.setdefault(str(r.chrom), []).append(
            (r.gene, int(r.pos_min) - PAD, int(r.pos_max) + PAD))
    chroms = set(windows)

    recs = []
    with gzip.open(paths.CLINVAR_VCF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom = line[:line.find("\t")]
            if chrom not in chroms:
                continue
            p = line.rstrip("\n").split("\t")
            pos = int(p[1])
            gene = next((gn for gn, lo, hi in windows[chrom] if lo <= pos <= hi), None)
            if gene is None:
                continue
            info = parse_info(p[7])
            recs.append(dict(gene=gene, chrom=chrom, pos=pos, ref=p[3], alt=p[4],
                             clinvar_id=p[2],
                             CLNSIG=info.get("CLNSIG", ""),
                             CLNREVSTAT=info.get("CLNREVSTAT", ""),
                             CLNVC=info.get("CLNVC", ""),
                             GENEINFO=info.get("GENEINFO", ""),
                             CLNDN=info.get("CLNDN", "")[:80]))
    cv = pd.DataFrame(recs)
    cv.to_csv(DERIV / "clinvar_raw_regions.csv", index=False)
    return cv


# 3. Add gold stars and the truth bucket
def label_clinvar(cv: pd.DataFrame) -> pd.DataFrame:
    cv = cv.copy()
    cv["gold_stars"] = cv.CLNREVSTAT.map(STARS).fillna(0).astype(int)
    cv["clin_label"] = cv.CLNSIG.map(clnsig_bucket)
    cv["chrom"] = cv.chrom.astype(str)
    cv.to_csv(DERIV / "clinvar_labeled_regions.csv", index=False)
    return cv


# 4. Load the assayed variants + their three-class call
def load_assay() -> pd.DataFrame:
    frames = []
    for g in GENES:
        d = pd.read_csv(OUT / f"{g}_annotated.csv",
                        usecols=["chrom", "pos", "ref", "alt", "route_class",
                                 "coarse_consequence", "protein_visible", "dist_bin", "score"])
        d["gene"] = g
        frames.append(d)
    assay = pd.concat(frames, ignore_index=True)
    assay["chrom"] = assay.chrom.astype(str)

    call = pd.read_csv(POOLED, usecols=["gene", "chrom", "pos", "ref", "alt", "call"])
    call["chrom"] = call.chrom.astype(str)
    assay = assay.merge(call, on=["gene", "chrom", "pos", "ref", "alt"], how="left")
    return assay


# 5. One ClinVar record per locus (best evidence), then join onto the assay
def join_clinvar(assay: pd.DataFrame, cv: pd.DataFrame) -> pd.DataFrame:
    # keep the strongest record at each locus: highest stars, then richest label
    lab_rank = {"P/LP": 0, "B/LB": 1, "conflicting": 2, "VUS": 3, "other": 4}
    cv = cv.assign(_lab=cv.clin_label.map(lab_rank))
    best = (cv.sort_values(["gold_stars", "_lab"], ascending=[False, True])
              .drop_duplicates(["chrom", "pos", "ref", "alt"], keep="first"))
    key = ["chrom", "pos", "ref", "alt"]
    j = assay.merge(best[key + ["clin_label", "gold_stars", "CLNSIG", "clinvar_id"]],
                    on=key, how="left")
    j.to_csv(DERIV / "assay_clinvar_joined.csv", index=False)
    return j


# 6. gnomAD combined allele frequency (benign cross-check)
def add_gnomad(j: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for g in GENES:
        payload = json.load(open(DERIV / f"gnomad_r4_{g}.json"))
        chrom = str(payload["chrom"])
        for v in payload["variants"]:
            ex = v.get("exome") or {}
            ge = v.get("genome") or {}
            ac = (ex.get("ac") or 0) + (ge.get("ac") or 0)
            an = (ex.get("an") or 0) + (ge.get("an") or 0)
            rows.append((chrom, v["pos"], v["ref"], v["alt"], ac, an, ac / an if an else 0.0))
    gdf = (pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt",
                                       "gnomad_ac", "gnomad_an", "gnomad_af"])
             .drop_duplicates(["chrom", "pos", "ref", "alt"]))
    gdf["chrom"] = gdf.chrom.astype(str)

    ctrl = j.merge(gdf, on=["chrom", "pos", "ref", "alt"], how="left")
    ctrl["in_gnomad"] = ctrl.gnomad_ac.notna() & (ctrl.gnomad_ac > 0)
    ctrl["gnomad_af"] = ctrl.gnomad_af.fillna(0.0)
    ctrl["gnomad_benign_proxy"] = ctrl.gnomad_af >= 1e-3     # combined AF, not popmax FAF
    return ctrl


# 7. Assign roles and freeze the fit/validate split
def lock_split(ctrl: pd.DataFrame) -> pd.DataFrame:
    is_ctrl = (ctrl.gold_stars >= 2) & ctrl.clin_label.isin(["P/LP", "B/LB"])
    ctrl["role"] = np.where(is_ctrl, "calibration", "apply")

    def fold(row):
        h = int(hashlib.md5(f"{row.chrom}:{row.pos}:{row.ref}:{row.alt}:{SEED}"
                            .encode()).hexdigest(), 16) % 1000
        return "fit" if h < 500 else "validate"

    ctrl["calib_fold"] = np.where(ctrl.role == "calibration", ctrl.apply(fold, axis=1), "")

    cols = ["chrom", "pos", "ref", "alt", "coarse_consequence", "route_class",
            "protein_visible", "score", "dist_bin", "gene", "call", "clin_label",
            "gold_stars", "CLNSIG", "clinvar_id", "gnomad_ac", "gnomad_an", "gnomad_af",
            "in_gnomad", "gnomad_benign_proxy", "role", "calib_fold"]
    ctrl[cols].to_csv(paths.CLINICAL_CONTROLS, index=False)
    return ctrl


# 8. Control-count tables (per gene, per gene x stratum)
def overlap_tables(ctrl: pd.DataFrame) -> None:
    def n(df, star, lab):
        return int(((df.gold_stars >= star) & (df.clin_label == lab)).sum())

    # per gene, at 1/2/3-star floors
    rows = []
    for g, d in ctrl.groupby("gene"):
        rows.append(dict(gene=g, assayed=len(d),
                         P_1star=n(d, 1, "P/LP"), B_1star=n(d, 1, "B/LB"),
                         P_2star=n(d, 2, "P/LP"), B_2star=n(d, 2, "B/LB"),
                         P_3star=n(d, 3, "P/LP"), B_3star=n(d, 3, "B/LB")))
    per_gene = pd.DataFrame(rows)
    pooled = per_gene.drop(columns="gene").sum()
    pooled["gene"] = "POOLED"
    per_gene = pd.concat([per_gene, pd.DataFrame([pooled])], ignore_index=True)
    per_gene.to_csv(DERIV / "overlap_per_gene.csv", index=False)

    # per gene x route_class
    a = ctrl.assign(
        P2=((ctrl.gold_stars >= 2) & (ctrl.clin_label == "P/LP")).astype(int),
        B2=((ctrl.gold_stars >= 2) & (ctrl.clin_label == "B/LB")).astype(int),
        P1=((ctrl.gold_stars >= 1) & (ctrl.clin_label == "P/LP")).astype(int),
        B1=((ctrl.gold_stars >= 1) & (ctrl.clin_label == "B/LB")).astype(int))
    st = (a.groupby(["gene", "route_class"])
            .agg(n_assayed=("pos", "size"),
                 P2=("P2", "sum"), B2=("B2", "sum"),
                 P1=("P1", "sum"), B1=("B1", "sum"))
            .reset_index())
    st.to_csv(DERIV / "overlap_per_stratum.csv", index=False)


# 9. Provenance + split manifest
def write_manifest(ctrl: pd.DataFrame) -> None:
    def counts(df):
        return {lab: int((df.clin_label == lab).sum()) for lab in ["P/LP", "B/LB"]}
    manifest = {
        "created_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "seed": SEED,
        "sources": {
            "clinvar_vcf": {"release": "2026-07-06", "file": "clinvar.vcf.gz",
                            "md5": "f78d25d49e17a070957a127e409f87b9", "build": "GRCh38"},
            "gnomad": {"dataset": "gnomad_r4", "via": "variants MCP region_variants"},
            "clinvar_mirror": {"release": "2026-06-06", "via": "variants MCP clinvar_variants"},
        },
        "join_key": "chrom,pos,ref,alt (GRCh38, VCF left-aligned)",
        "label_rule": {
            "P/LP": sorted(P_SET), "B/LB": sorted(B_SET),
            "primary_control_star_floor": 2, "sensitivity_tier_star_floor": 1,
            "gnomad_benign_proxy": "combined AF>=1e-3 (exome+genome; popmax FAF not in payload)"},
        "roles": {
            "calibration": "2*+ ClinVar P/LP or B/LB; fit/validate 50/50 by frozen md5 hash",
            "apply": "all other assayed variants (VUS/conflicting/no-ClinVar/low-star) "
                     "= reclassification targets, never inform thresholds"},
        "role_counts": {
            "calibration": int((ctrl.role == "calibration").sum()),
            "apply": int((ctrl.role == "apply").sum()),
            "calibration_fit": int((ctrl.calib_fold == "fit").sum()),
            "calibration_validate": int((ctrl.calib_fold == "validate").sum())},
        "controls_2star_per_gene": {
            g: counts(d[d.role == "calibration"]) for g, d in ctrl.groupby("gene")},
    }
    json.dump(manifest, open(paths.SPLIT_MANIFEST, "w"), indent=2)


def main() -> None:
    reg   = gene_regions()
    cv    = label_clinvar(parse_clinvar(reg))
    assay = load_assay()
    ctrl   = lock_split(add_gnomad(join_clinvar(assay, cv)))
    overlap_tables(ctrl)
    write_manifest(ctrl)
    print(f"built clinical_controls.csv: {len(ctrl)} variants | "
          f"calibration={int((ctrl.role=='calibration').sum())} "
          f"apply={int((ctrl.role=='apply').sum())}")


if __name__ == "__main__":
    main()
