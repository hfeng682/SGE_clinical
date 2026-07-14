"""

Layout:
    data/            inputs (screens, pooled calls, clinical controls, AG scores)
    data_processing/ raw inputs -> analysis-ready tables
    model/           analysis-ready tables -> variant evidence codes
    validation/      out-of-sample checks
    results/         everything the model and validation write
    manuscript/      figures/ and tables/ for the write-up
"""
from __future__ import annotations
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]          # hackathon/

# ---- data (inputs) --------------------------------------------------------
DATA              = REPO / "data"
RAW_SCREENS       = DATA / "raw_screens"             # MaveDB score files
ANNOTATED_SCREENS = DATA / "annotated_screens"       # VEP-annotated screens
POOLED            = DATA / "pooled_labeled.csv"      # shared 3-class call per variant
CLINICAL_CONTROLS = DATA / "clinical_controls.csv"   # ClinVar/gnomAD-labelled controls
SPLIT_MANIFEST    = DATA / "split_manifest.json"     # frozen calibration/apply split
CLINVAR_GNOMAD    = DATA / "clinvar_gnomad"          # control sources + derivation tables
CLINVAR_VCF       = CLINVAR_GNOMAD / "clinvar.vcf.gz"
AG_SCORES         = DATA / "alphagenome_scores.csv"  # AlphaGenome per-variant scores
AG_SHARDS         = DATA / "alphagenome_shards"

# ---- results (outputs) ----------------------------------------------------
RESULTS             = REPO / "results"
RESULTS_BOUNDS      = RESULTS / "bounds"             # error-bound engine outputs
RESULTS_STRENGTH    = RESULTS / "strength"           # OddsPath strength outputs
RESULTS_CODES       = RESULTS / "evidence_codes"     # final per-variant codes
RESULTS_RELIABILITY = RESULTS / "reliability_map"    # AlphaGenome reliability map
RESULTS_VALIDATION  = RESULTS / "validation"

# ---- manuscript -----------------------------------------------------------
MANUSCRIPT = REPO / "manuscript"
FIGURES    = MANUSCRIPT / "figures"
TABLES     = MANUSCRIPT / "tables"

# ---- the five in-session screens ------------------------------------------
GENES = ["BARD1", "BRCA1", "PALB2", "RAD51D", "VHL"]

for _d in (RESULTS, RESULTS_BOUNDS, RESULTS_STRENGTH, RESULTS_CODES,
           RESULTS_RELIABILITY, RESULTS_VALIDATION, FIGURES, TABLES):
    _d.mkdir(parents=True, exist_ok=True)
