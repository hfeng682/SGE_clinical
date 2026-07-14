"""
tissue.py -- Stage 4: the tissue-transfer discount.

A loss that acts through the tissue-INVARIANT protein route (a broken protein
product) transfers to a patient's disease tissue reliably; a loss that acts
through the tissue-VARIABLE RNA route (splicing / abundance) does not, unless
it is corroborated in a disease-relevant tissue. The method encodes this by
dropping an uncorroborated RNA-route loss exactly ONE ACMG tier
(Strong -> Moderate -> Supporting -> uncertain); a protein-route loss is never
discounted.

The one-tier size is validated in validation/tissue_discount.py.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

# BRCA1/VHL carry in-hand cross-conditions (independent cell-types / timepoints).
# The frozen rna_drop5 flag is driven by:  BARD1/PALB2/RAD51D -> rna_score,
# BRCA1 -> score_rna (combined),  VHL -> rna_score_d20.  The OTHER conditions
# are the independent corroboration sources used here.
CROSS_CONDITION_COLS = {
    "BRCA1": ["score_rna_rep1", "score_rna_rep2"],   # two cell-type replicates
    "VHL":   ["rna_score_d6"],                         # the other timepoint (d20 is primary)
}
# preliminary FPR-5% RNA-drop thresholds (the cut that DEFINES the frozen flag)
RNA_THR5 = {"BARD1": -0.684, "PALB2": -0.564, "RAD51D": -1.016,
            "BRCA1": -0.667, "VHL": -0.597}


def tissue_corroboration(r: pd.Series) -> tuple[bool, str]:
    """Is this RNA-route loss corroborated in disease-relevant tissue?

    Three admissible sources of corroboration (any one suffices):
      A. in-hand cross-condition: the mRNA drop reproduces in an INDEPENDENT
         condition of the same screen (BRCA1 both cell-type replicates; VHL the
         other timepoint) -> the effect is condition-invariant.
      B. AG-gated splice: for an unmeasurable variant inside the AlphaGenome
         junction-proximal reliability envelope, AG predicts a splice disruption
         whose disease-tissue expression is concordant and near tissue-invariant
         (transfer_verdict == 'splice_disruptor_expr_concordant').
      C. measured & already tissue-invariant is NOT auto-granted -- a single
         measured drop in the assay's own cell line is the observation being
         transferred, not corroboration of transfer. It needs A or B.

    Returns (corroborated, corroboration_source).
    """
    gene = r["gene"]
    thr = RNA_THR5[gene]
    # --- source A: cross-condition reproduction (BRCA1 / VHL) ---
    cols = CROSS_CONDITION_COLS.get(gene, [])
    if cols and bool(r.get("rna_measured", False)):
        vals = [r.get(c, np.nan) for c in cols]
        vals = [v for v in vals if pd.notna(v)]
        if vals and all(v <= thr for v in vals):
            return True, f"cross_condition({'+'.join(cols)})"
    # --- source B: AG-gated splice disruptor, concordant disease tissue ---
    if r.get("transfer_verdict") == "splice_disruptor_expr_concordant":
        return True, "ag_splice_gated(disease_concordant)"
    return False, "none"
