"""
mechanism.py -- assign each variant a mechanism and its ACMG/AMP strength.

Three stages, no new statistics -- every tier is read from the frozen OddsPath
strength tables (model/build_strength.py):

  STAGE 1  mechanism         : each variant -> a mechanism (protein / RNA / ...)
  STAGE 2  strength ladder   : mechanism + gene x class -> an error-bounded ACMG tier
  STAGE 3  SVI code          : mechanism + direction    -> a ClinGen SVI code

The tissue-transfer discount (Stage 4) lives in tissue.py; final assembly in
build_evidence_codes.py.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # repo root (utils)
from utils import paths

# ACMG strength ladder, ordered. "uncertain" is the floor (no transferable code).
TIER_ORDER = ["uncertain", "Supporting", "Moderate", "Strong"]
TIER_RANK = {t: i for i, t in enumerate(TIER_ORDER)}


def downweight(tier: str, n: int = 1) -> str:
    """Move DOWN the ACMG ladder by n rungs; floor at 'uncertain'."""
    return TIER_ORDER[max(0, TIER_RANK.get(tier, 0) - n)]


# Mechanistic class partition.
# NB: the raw `protein_visible` column means "sits in the CDS / has a protein
# coordinate" -- it is TRUE for synonymous variants.  The mechanistic
# protein-altering split is `protein_altering` (== route_class in PROTEIN_ALT).
PROTEIN_ALT = {"missense", "inframe/indel", "nonsense", "start/stop_lost"}
PROTEIN_BLIND_EXONIC = {"synonymous", "splice_region"}   # measurable -> RNA route
PROTEIN_BLIND_INTRON = {"intronic"}                        # logic-assigned RNA route
UTR = {"3UTR", "5UTR"}


# Stage 1: mechanism assignment
def route_mechanism(r: pd.Series) -> str:
    """Assign one mechanism to a variant from its frozen call + flags.

    Mechanisms
    ----------
    uncertain     : the survival call is Uncertain (guaranteed abstention).
    benign        : the survival call is Normal (benign direction; tissue-invariant,
                    never discounted -- a normal transcript transfers everywhere).
    protein       : protein-altering LoF with no measured mRNA drop
                    (broken protein, intact transcript -> tissue-invariant).
    dual          : protein-altering LoF that also drops mRNA
                    (protein + RNA components compete -- the 501 cohort).
                    NB nonsense mRNA drop is NMD, the expected consequence of a
                    premature stop, so nonsense stays 'protein', not 'dual'.
    rna_measured  : protein-blind exonic LoF with a measured mRNA drop.
    rna_nodrop    : protein-blind exonic LoF, measurable, but NO mRNA drop
                    (survival calls loss; the RNA score does not corroborate).
    rna_logic     : protein-blind non-measurable LoF (intronic/deep-splice) --
                    transcript-mediated by construction, survival supplies the LoF.
    utr           : UTR LoF (judged against its own abundance baseline).
    """
    call = r["call"]
    if call == "Uncertain":
        return "uncertain"
    if call == "Normal":
        return "benign"
    # ---- call == 'LoF' ----
    rc = r["route_class"]
    if rc in UTR:
        return "utr"
    drop = bool(r["rna_drop5"])
    meas = bool(r["rna_measured"])
    if rc in PROTEIN_ALT:
        # nonsense drop is NMD (expected) -> protein mechanism, not dual
        if rc == "nonsense":
            return "protein"
        return "dual" if (meas and drop) else "protein"
    if rc in PROTEIN_BLIND_EXONIC:
        if meas and drop:
            return "rna_measured"
        if meas and not drop:
            return "rna_nodrop"
        return "rna_logic"           # protein-blind exonic but not RNA-measurable
    if rc in PROTEIN_BLIND_INTRON:
        return "rna_logic"
    return "uncertain"               # unreachable given the partition


# mechanism -> which ACMG direction it speaks to.
# rna_nodrop IS a LoF call (survival called loss) on a protein-blind variant:
# the protein is unchanged, so the loss is transcript-mediated BY CONSTRUCTION
# even though the abundance score showed no drop -- it is the RNA mechanism with
# the weakest transfer evidence, and is (almost always) uncorroborated below.
ROUTE_IS_LOF = {
    "protein": True, "dual": True, "rna_measured": True,
    "rna_logic": True, "utr": True, "rna_nodrop": True,
    "benign": False, "uncertain": False,
}
# whether the mechanism is a transcript/RNA route (subject to the tissue discount)
ROUTE_IS_RNA = {
    "rna_measured": True, "rna_logic": True, "utr": True, "rna_nodrop": True,
    "dual": True,          # dual: functional tier kept, transfer discounted
    "protein": False, "benign": False, "uncertain": False,
}


# Stage 3: SVI code
def svi_code(route: str, direction: str) -> str:
    """ClinGen SVI splicing-code partition.

    Protein mechanism  -> PS3 (LoF) / BS3 (normal).
    RNA mechanism      -> PVS1_RNA (LoF via transcript) / BP7 (no transcript impact).
    """
    if route in ("protein",):
        return "PS3" if direction == "LoF" else "BS3"
    if route in ("dual",):
        # dual carries both a protein call and an RNA component; the functional
        # code is PS3 (protein evidence present), annotated with the RNA flag.
        return "PS3"
    if route in ("rna_measured", "rna_logic", "rna_nodrop", "utr"):
        return "PVS1_RNA" if direction == "LoF" else "BP7"
    if route == "benign":
        return "BS3"
    return "none"


# Stage 2: strength ladder (error-bounded, per gene x class x direction)
MIN_N = 10   # frozen per-side control floor (matches build_strength MIN_N)


def load_strength_sources():
    """Assemble the error-bounded strength lookup from the OddsPath tables.

    Returns two lookups:
      per_class[(gene, route_class, direction)] ->
          dict(tier, oddspath_lo, abstain, abstain_reason, n_P, n_B, evaluable)
        [ALL groups incl. abstained; `evaluable` = both sides >= MIN_N, i.e. the
         class was genuinely testable in this gene]
      per_gene[(gene, direction)]               -> (tier, oddspath_lo)
    """
    pc = pd.read_csv(paths.RESULTS_STRENGTH / "oddspath_calls_per_class.csv")
    per_class = {}
    for r in pc.itertuples():
        evaluable = (r.n_P >= MIN_N) and (r.n_B >= MIN_N)
        per_class[(r.gene, r.route_class, r.direction)] = dict(
            tier=r.tier, oddspath_lo=r.oddspath_lo,
            abstain=bool(r.abstain),
            abstain_reason=("" if pd.isna(r.abstain_reason) else r.abstain_reason),
            n_P=int(r.n_P), n_B=int(r.n_B), evaluable=evaluable)
    pg = pd.read_csv(paths.RESULTS_STRENGTH / "oddspath_per_gene.csv")
    per_gene = {
        (r.gene, r.direction): (r.tier, r.oddspath_lo)
        for r in pg[pg.abstain == False].itertuples()
    }
    return per_class, per_gene


# the ACMG call direction (LoF / Normal) reads the corresponding OddsPath ladder:
#   a LoF call is scored on the PS3 (pathogenic) ladder,
#   a Normal call is scored on the BS3 (benign) ladder.
DIRECTION_TO_ODDSPATH = {"LoF": "PS3", "Normal": "BS3"}


def assign_strength(gene, route_class, direction, route,
                    per_class, per_gene):
    """Error-bounded ACMG tier for a LoF/normal call, with its provenance.

    `direction` is the ACMG call direction ('LoF' or 'Normal'); it is mapped to
    the OddsPath side ('PS3' / 'BS3') the strength tables are keyed on.

    Hierarchy (respecting the frozen single-gene / single-class discipline):

      1. CLASS-MATCHED, genuinely evaluated (both control sides >= MIN_N):
         use the class result directly. If that group was error-bounded, take
         its tier; if it was EVALUATED but abstained (e.g. VHL splice BS3, bound
         below Supporting), respect the abstention -> uncertain. This is the only
         path for splice-region (all 5 genes) and BRCA1 missense.

      2. NOT class-evaluable (no both-sided class controls) -> whole-gene fallback,
         but only where it is defensible:
           * PROTEIN-ALTERING class (missense/nonsense/inframe/start-stop): take
             the whole-gene OddsPath (the whole-gene number is a same-kind protein
             proxy). Provenance 'gene_oddspath' carries the per-class caveat.
           * PROTEIN-BLIND class (synonymous/intronic/UTR):
               - BENIGN (BS3) direction: whole-gene BS3 is supported by exactly
                 these benign controls -> use it ('gene_oddspath_benign').
               - PATHOGENIC (PS3) direction: NO class-matched pathogenic control
                 exists, and mixing protein classes into a protein-blind LoF call
                 is not defensible -> abstain to uncertain clinical tier. The call
                 still carries its survival evidence and its mechanism; it simply
                 gets no class clinical PS3.

    Returns (tier, oddspath_lo, strength_source).
    """
    op_dir = DIRECTION_TO_ODDSPATH.get(direction)
    if op_dir is None:
        return "uncertain", np.nan, "abstain_no_control"

    row = per_class.get((gene, route_class, op_dir))
    if row is not None and row["evaluable"]:
        # class was genuinely testable in this gene
        if not row["abstain"]:
            return row["tier"], row["oddspath_lo"], "class_oddspath"
        return "uncertain", np.nan, "class_abstained:" + row["abstain_reason"]

    # not class-evaluable -> whole-gene fallback, gated on mechanism class
    if route_class in PROTEIN_ALT:
        hit = per_gene.get((gene, op_dir))
        if hit is not None:
            return hit[0], hit[1], "gene_oddspath"
        return "uncertain", np.nan, "abstain_no_gene_fallback"
    # protein-blind class
    if op_dir == "BS3":
        hit = per_gene.get((gene, op_dir))
        if hit is not None:
            return hit[0], hit[1], "gene_oddspath_benign"
        return "uncertain", np.nan, "abstain_no_gene_fallback"
    # protein-blind, pathogenic direction: no defensible class clinical strength
    return "uncertain", np.nan, "abstain_protein_blind_no_class_pathogenic_control"
