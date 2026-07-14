"""
strength.py -- clinical OddsPath strength with a finite-sample error bound.

Turns the frozen three-class assay call into an ACMG/AMP functional-evidence
strength, using the clinical controls (ClinVar/gnomAD) and the Clopper-Pearson error-bound
engine UNCHANGED. Nothing here fits a threshold; the survival-score call is
already frozen (investigator cuts + frozen RNA-drop cuts), so this module
only maps an already-made call onto the OddsPath ladder.

KEY IDENTITY (why the binary-error engine applies without modification)
------------------------------------------------------------------------------
With P1 taken as the empirical pathogenic proportion among a group's controls
(the project's convention, step0_analysis.py), the ClinGen/Tavtigian OddsPath
collapses to a likelihood ratio on two DISJOINT, INDEPENDENT control samples:

    OddsPath_PS3 = [P2 (1-P1)] / [(1-P2) P1]
                 = (a / n_P) / (b / n_B)
                 = TPR / FPR

    a = pathogenic controls called LoF,  n_P = pathogenic controls
    b = benign     controls called LoF,  n_B = benign controls

TPR is a Bernoulli proportion on the n_P pathogenic controls; FPR is a Bernoulli
proportion on the n_B benign controls. They are independent because the control
sets are disjoint. Each is exactly the {0,1} binary-call-error quantity the error-bound stage
earns a code with Clopper-Pearson on the position-clustered effective sample. The
conservative (lower) OddsPath is therefore

    OddsPath_PS3_lo = TPR_lo / FPR_hi

with TPR_lo a lower CP bound and FPR_hi an upper CP bound, each at confidence
1 - delta/2 so the pair holds jointly at >= 1 - delta (independent samples;
Bonferroni). The tier is read off this conservative bound.

BS3 (benign) direction is symmetric on the NORMAL call:

    OddsPath_BS3 = (benign called Normal / n_B) / (path called Normal / n_P)
                 = TNR / FNR ,   conservative  OddsPath_BS3_lo = TNR_lo / FNR_hi

Both directions read on the SAME 2.1 / 4.3 / 18.7 ladder (larger = stronger).
The conventional ClinGen "benign OddsPath" is the reciprocal 1/OddsPath_BS3.

FROZEN RULES (inherited from Steps 0-1, enforced here)
------------------------------------------------------
  * Every group is one gene x one variant class. Classes are never merged and
    genes are never pooled. A per-class OddsPath needs both pathogenic and
    benign controls in that gene x class; single-sided classes (nonsense = path
    only; synonymous/intronic/UTR = benign only) get NO OddsPath here and are
    covered by the error-bound stage instead.
  * Only the calibration controls (role == 'calibration'; 2-star ClinVar P/LP or
    B/LB) ever enter. Apply variants never touch a number in this module.
  * Abstain (return "uncertain", no code) when either side has n < MIN_N or the
    conservative bound does not clear the Supporting rung.

Self-contained beyond bounds.py (numpy + scipy + pandas).
"""

from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import beta as _beta

# The shared error-bound engine lives beside this module.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bounds import ub_cp, _cluster_cp  # error-bound engine, reused unchanged

# ACMG/AMP OddsPath ladder (Tavtigian 2018 prior 0.10; ClinGen SVI / Brnich).  #
# PS3 pathogenic thresholds; BS3 read on the same scale via the reciprocal.    #
ODDS_SUPPORTING, ODDS_MODERATE, ODDS_STRONG = 2.1, 4.3, 18.7
DELTA = 0.05          # overall per-group-per-direction target: 1 - DELTA = 0.95
MIN_N = 10            # frozen abstention floor (per side)


def tier_from_oddspath(op: float) -> str:
    """Strength rung for a conservative OddsPath on the shared ladder."""
    if not np.isfinite(op):
        return "n/a"
    if op >= ODDS_STRONG:
        return "Strong"
    if op >= ODDS_MODERATE:
        return "Moderate"
    if op >= ODDS_SUPPORTING:
        return "Supporting"
    return "below"


# Lower-tail Clopper-Pearson (the only addition to the error-bound engine).         #
# ub_cp already gives the upper tail; lb_cp is its exact mirror.               #

def lb_cp(k, n, delta: float = 0.05) -> float:
    """One-sided LOWER Clopper-Pearson bound: BetaInv(delta; k, n-k+1).

    Exact for i.i.d. Bernoulli; accepts non-integer effective (k, n) for the
    position-clustered plug-in, exactly as ub_cp does for the upper tail.
    """
    if n <= 0:
        return 0.0
    if k <= 0:
        return 0.0
    if k >= n:
        # all-success: lower bound is BetaInv(delta; n, 1)
        return float(_beta.ppf(delta, n, 1.0))
    return float(_beta.ppf(delta, k, n - k + 1.0))


def _cluster_counts(err: np.ndarray, pos: np.ndarray) -> tuple[int, float]:
    """(n_clusters, k_clusters) from position-clustering: one independent draw
    per genomic position, per-position mean error. Mirrors bounds._cluster_cp
    but returns the effective counts so both CP tails can be taken."""
    per = (pd.DataFrame({"pos": np.asarray(pos), "e": np.asarray(err, float)})
           .groupby("pos")["e"].mean())
    return int(len(per)), float(per.sum())


def cp_rate_two_sided(err: np.ndarray, pos: np.ndarray,
                      delta_each: float) -> dict:
    """Position-clustered CP for a binary indicator, both tails, at 1-delta_each.

    Returns raw n/k, cluster n/k, point rate, and clustered lower & upper CP.
    Consistency check: the upper equals bounds._cluster_cp exactly.
    """
    err = np.asarray(err, float)
    pos = np.asarray(pos)
    n = int(len(err))
    if n == 0:
        return dict(n=0, k=0.0, n_cl=0, k_cl=0.0, point=np.nan,
                    lo=0.0, hi=1.0)
    k = float(err.sum())
    n_cl, k_cl = _cluster_counts(err, pos)
    lo = lb_cp(k_cl, n_cl, delta_each)
    hi = ub_cp(k_cl, n_cl, delta_each)
    return dict(n=n, k=k, n_cl=n_cl, k_cl=k_cl, point=k / n, lo=lo, hi=hi)


# OddsPath for one gene x class group, one direction, with a error-bounded bound.

def oddspath_group(sub: pd.DataFrame, direction: str, delta: float = DELTA,
                   min_n: int = MIN_N, pos_col: str = "pos") -> dict:
    """Error-bounded OddsPath for one single-gene, single-class control group.

    sub        : calibration controls for one gene x one class, with columns
                 clin_label in {P/LP, B/LB}, `call` in {LoF, Uncertain, Normal},
                 and a genomic-position column `pos_col`.
    direction  : "PS3" (pathogenic, from the LoF call) or
                 "BS3" (benign, from the Normal call).

    Both directions need both pathogenic and benign controls; the strength is
    read off the conservative end of a position-clustered Clopper-Pearson bound,
    split delta/2 per independent side. Returns a flat dict (one table row).
    """
    delta_each = delta / 2.0
    P = sub[sub.clin_label == "P/LP"]
    B = sub[sub.clin_label == "B/LB"]
    n_P, n_B = len(P), len(B)

    rec = dict(direction=direction, n_P=n_P, n_B=n_B)

    if direction == "PS3":
        # numerator lives on the pathogenic side (TPR, want lower);
        # denominator on the benign side (FPR, want upper).
        num_err, num_pos = (P["call"] == "LoF").values, P[pos_col].values
        den_err, den_pos = (B["call"] == "LoF").values, B[pos_col].values
        num_name, den_name = "TPR", "FPR"
    elif direction == "BS3":
        # numerator on the benign side (TNR, want lower);
        # denominator on the pathogenic side (FNR, want upper).
        num_err, num_pos = (B["call"] == "Normal").values, B[pos_col].values
        den_err, den_pos = (P["call"] == "Normal").values, P[pos_col].values
        num_name, den_name = "TNR", "FNR"
    else:
        raise ValueError("direction must be 'PS3' or 'BS3'")

    num = cp_rate_two_sided(num_err, num_pos, delta_each)
    den = cp_rate_two_sided(den_err, den_pos, delta_each)

    rec.update({
        f"{num_name}_point": num["point"], f"{num_name}_lo": num["lo"],
        f"{den_name}_point": den["point"], f"{den_name}_hi": den["hi"],
        "n_cl_num": num["n_cl"], "n_cl_den": den["n_cl"],
    })

    # point OddsPath (report inf when the denominator's point rate is 0)
    if den["point"] and den["point"] > 0:
        rec["oddspath_point"] = num["point"] / den["point"]
    else:
        rec["oddspath_point"] = np.inf if num["point"] and num["point"] > 0 else np.nan

    # abstention on control counts (both sides needed)
    if n_P < min_n or n_B < min_n:
        rec.update(oddspath_lo=np.nan, tier="uncertain", abstain=True,
                   abstain_reason=("n_P<%d" % min_n if n_P < min_n
                                   else "n_B<%d" % min_n))
        return rec

    # conservative OddsPath = lower_numerator / upper_denominator (finite: den hi > 0)
    op_lo = num["lo"] / den["hi"] if den["hi"] > 0 else np.inf
    rec["oddspath_lo"] = op_lo
    t = tier_from_oddspath(op_lo)
    rec["tier"] = t if t != "below" else "uncertain"
    rec["abstain"] = (t == "below")
    rec["abstain_reason"] = "bound_below_supporting" if t == "below" else ""
    return rec


def oddspath_grid(controls: pd.DataFrame, gene_col: str = "gene",
                  class_col: str = "route_class", delta: float = DELTA,
                  min_n: int = MIN_N, pos_col: str = "pos") -> pd.DataFrame:
    """Error-bound the OddsPath for every gene x class group, both directions.

    controls : the calibration set (role == 'calibration'), already restricted
               to whatever fold the caller wants (fit / validate / full).
    Returns one row per (gene, class, direction).
    """
    rows = []
    for (g, cl), sub in controls.groupby([gene_col, class_col], observed=True):
        for direction in ("PS3", "BS3"):
            rec = dict(gene=g, route_class=cl)
            rec.update(oddspath_group(sub, direction, delta, min_n, pos_col))
            rows.append(rec)
    lead = ["gene", "route_class", "direction", "n_P", "n_B",
            "oddspath_point", "oddspath_lo", "tier", "abstain", "abstain_reason"]
    df = pd.DataFrame(rows)
    rest = [c for c in df.columns if c not in lead]
    return df[lead + rest].sort_values(["gene", "route_class", "direction"]
                                       ).reset_index(drop=True)
