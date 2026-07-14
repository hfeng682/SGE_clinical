"""
analysis_common.py — shared definitions for the reliability map.

Everything downstream of `alphagenome_per_variant.csv` imports from here so the
reliability map, tissue-bridge, GTEx-audit and figure scripts share one frozen set of
definitions. Mirrors `ag_common.py` on the scoring side.

Main story (start to end):
  scoring (ag_common/run_shard/combine_shards) -> alphagenome_per_variant.csv
     |
  reliability.py         reliability of AG vs the MEASURED mRNA-drop event, by distance
  tissue_bridge.py AG tissue evidence for UNMEASURABLE variants, gated to the
                   distance zone the reliability map proved reliable
  gtex_audit.py    AG tissue tracks ARE GTEx v8 -> do not double-count vs sQTL
  figures.py       the one figure that carries the finding
"""
import os
import numpy as np
import pandas as pd

# ---- paths -----------------------------------------------------------------
import sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
_sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))   # repo root
from utils import paths
DEFAULT_INPUT = str(paths.AG_SCORES)                          # data/alphagenome_scores.csv
RESULTS_DIR = str(paths.RESULTS_RELIABILITY)                  # results/reliability_map/

# ---- ground truth (frozen) -------------------------------------------------
# The SGE assay measures mRNA abundance for EXONIC variants only. Those are the
# reliability map ground truth; the unmeasurable intronic variants are the tissue-bridge
# targets, disjoint by construction.
MEASURED_FLAG = "rna_measured"          # True == exonic, measurable
EVENT = "rna_drop5"                     # primary event: measured mRNA-abundance drop (5% threshold)

# ---- AlphaGenome predictors ------------------------------------------------
SPLICE_MAG_COLS = ["splice_sites_maxabs", "splice_usage_maxabs", "splice_junctions_maxabs"]
EXPR_COL = "expr_lfc_maxabs"            # gene-level expression LFC magnitude
EXPR_SIGNED_COL = "expr_lfc_signed_at_maxabs"

def ag_splice_mag(df):
    """Per-variant AG splice magnitude = max of the three splice scorers.
    (Raw scales differ; junction-dominated. Each is also evaluated alone in reliability.py.)"""
    return df[SPLICE_MAG_COLS].max(axis=1)

# ---- distance strata -------------------------------------------------------
DIST_BIN = "dist_bin"
DIST = "dist_to_junction"
# fine bins within the exon core, on distance to junction
FINE_EDGES = [3, 10, 25, 50, 100, np.inf]
FINE_LABELS = ["4-10", "11-25", "26-50", "51-100", ">100"]

def fine_bin(dist):
    """Label a distance (nt) into the within-exon-core fine bins; None if <=3."""
    d = np.asarray(dist, dtype=float)
    out = np.full(d.shape, None, dtype=object)
    for lo, hi, lab in zip([3, 10, 25, 50, 100], FINE_EDGES[1:], FINE_LABELS):
        out[(d > lo) & (d <= hi)] = lab
    return out

# ---- reliability rule (frozen) ---------------------------------------------
# A stratum is reliable when the best of the splice scores has bootstrap
# AUROC lower-CI > 0.70. The reliability map shows this holds ONLY at the junction edge
# (<=3 nt); reliability collapses to chance mid-exon.
RELIABLE_CI_LO = 0.70

# The reliable zone as a set of dist_bin labels (junction-proximal, <=~8 nt of a
# junction). By distance symmetry these are the strata a tissue-bridge inference
# may inherit reliability from. Exact strings must match the CSV.
RELIABLE_DIST_BINS = {
    "exonic_splice_region(\u22643)",     # exonic_splice_region(<=3)  -- the reliable reliability map stratum
    "splice_site(\u00b11-2)",            # splice_site(+/-1-2)         -- intronic, symmetric
    "splice_region_intronic(3-8)",       # intronic, symmetric
}

# ---- operating points ------------------------------------------------------
SPEC_TARGETS = [0.90, 0.95]             # specificity levels for missed-demotion counts

# ---- bootstrap -------------------------------------------------------------
BOOTSTRAP_N = 2000
SEED = 0
MIN_N, MIN_EVENTS = 20, 5               # power filter: skip thinner strata

# ---- disease-relevant tissues (for the tissue bridge) ----------------------
# HBOC genes -> breast/ovary/fallopian; VHL -> kidney/adrenal/CNS.
DISEASE_TISSUES = {
    "BRCA1":  ["Breast_Mammary_Tissue", "Ovary", "Fallopian_Tube"],
    "BARD1":  ["Breast_Mammary_Tissue", "Ovary", "Fallopian_Tube"],
    "PALB2":  ["Breast_Mammary_Tissue", "Ovary", "Fallopian_Tube"],
    "RAD51D": ["Breast_Mammary_Tissue", "Ovary", "Fallopian_Tube"],
    "VHL":    ["Kidney_Cortex", "Kidney_Medulla", "Adrenal_Gland", "Brain_Cerebellum", "Brain_Cortex"],
}

# ---- GTEx v8 tissue panel
GTEX_V8_TISSUES = [
    "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta",
    "Artery_Coronary", "Artery_Tibial", "Bladder", "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue",
    "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Cervix_Ectocervix",
    "Cervix_Endocervix", "Colon_Sigmoid", "Colon_Transverse",
    "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis",
    "Fallopian_Tube", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex",
    "Kidney_Medulla", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal",
    "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate",
    "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg",
    "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus",
    "Vagina", "Whole_Blood",
]

# ---- loading ---------------------------------------------------------------
def load_table(path=DEFAULT_INPUT):
    df = pd.read_csv(path, low_memory=False)
    assert (df["ag_status"] == "ok").all(), "some variants lack AG scores"
    return df

def measurable(df):
    """The reliability map ground-truth frame: measurable exonic variants."""
    return df[df[MEASURED_FLAG] == True].copy()

def unmeasurable(df):
    """The tissue-bridge target frame: intronic/unmeasurable variants."""
    return df[df[MEASURED_FLAG] == False].copy()

# ---- discrimination metrics (no sklearn dependency) ------------------------
def auroc(y, s):
    """AUROC via the Mann-Whitney rank statistic (average-rank ties)."""
    y = np.asarray(y).astype(bool)
    s = np.asarray(s, dtype=float)
    ok = ~np.isnan(s)
    y, s = y[ok], s[ok]
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    sr = s[order]
    i = 0
    while i < len(sr):
        j = i
        while j + 1 < len(sr) and sr[j + 1] == sr[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0  # 1-based average rank
        i = j + 1
    return (ranks[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)

def auprc(y, s):
    """Area under the precision-recall curve (trapezoidal; replicates sklearn PR+auc)."""
    y = np.asarray(y).astype(bool)
    s = np.asarray(s, dtype=float)
    ok = ~np.isnan(s)
    y, s = y[ok], s[ok]
    P = int(y.sum())
    if P == 0:
        return np.nan
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    tp = np.cumsum(y)
    fp = np.cumsum(~y)
    prec = tp / (tp + fp)
    rec = tp / P
    rec = np.concatenate([[0.0], rec])
    prec = np.concatenate([[1.0], prec])
    _trap = getattr(np, "trapezoid", None) or np.trapz  # numpy>=2 renamed trapz
    return float(_trap(prec, rec))

def bootstrap_ci(y, s, n=BOOTSTRAP_N, seed=SEED, alpha=0.05):
    """Percentile bootstrap 95% CI for AUROC."""
    y = np.asarray(y).astype(bool)
    s = np.asarray(s, dtype=float)
    ok = ~np.isnan(s)
    y, s = y[ok], s[ok]
    if y.sum() < 1 or (~y).sum() < 1:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    N = len(y)
    vals = np.empty(n)
    for b in range(n):
        idx = rng.integers(0, N, N)
        vals[b] = auroc(y[idx], s[idx])
    vals = vals[~np.isnan(vals)]
    return (float(np.percentile(vals, 100 * alpha / 2)),
            float(np.percentile(vals, 100 * (1 - alpha / 2))))

def threshold_at_specificity(neg_scores, spec):
    """Score threshold achieving `spec` specificity on the negative (non-event) scores."""
    neg = np.asarray(neg_scores, dtype=float)
    neg = neg[~np.isnan(neg)]
    return float(np.quantile(neg, spec))
