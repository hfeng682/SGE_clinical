"""
Finite-sample error bounds for the survival-score call.

Design rule: MATCH THE BOUND TO THE QUANTITY.

  * Binary call errors (false-LoF on synonymous negatives; false-Normal on
    nonsense positives) are exactly Bernoulli 0/1, so the error-bounded ceiling is
    the exact-binomial CLOPPER-PEARSON bound, applied to the correlation-aware
    effective sample (position-clustering is the primary widening; block
    bootstrap and a design-effect/VIF plug-in are the stability check). Raw-count
    CP is never reported as the error-bounded ceiling, because nearby variants are
    not independent.

  * The data-chosen RNA-drop threshold is an order statistic of the controls, so
    its error-bounded ceiling is the WILKS distribution-free order-statistic
    tolerance bound.

Self-contained (numpy + scipy only) so the same code is reused by the strength
module.
"""

from __future__ import annotations
import math
import numpy as np
import pandas as pd
from scipy.stats import beta as _beta

# Clopper-Pearson exact binomial — LEAD for binary call errors                #

def ub_cp(k, n, delta: float = 0.05) -> float:
    """One-sided upper Clopper-Pearson bound: BetaInv(1 - delta; k+1, n-k).

    Exact for i.i.d. Bernoulli.  Accepts non-integer (k, n) for the
    effective-sample-size plug-in used by the correlation widening.
    """
    if n <= 0:
        return 1.0
    if k >= n:
        return 1.0
    return float(_beta.ppf(1.0 - delta, k + 1.0, n - k))


# Wilks order-statistic tolerance bound — LEAD for the RNA-drop threshold      

# The threshold is thr = X_(k), the k-th order statistic of n control draws, so
# the true exceedance prob F(thr) ~ Beta(k, n-k+1) EXACTLY for any continuous F:
#     FPR <= BetaInv(1 - delta ; k, n-k+1)   w.p. >= 1 - delta.
# Exact, distribution-free, and correctly charges for the threshold having been
# chosen from the same controls.

def wilks_tolerance_upper(k: int, n: int, delta: float = 0.05) -> float:
    """Exact distribution-free upper bound on the exceedance prob of the k-th
    order statistic of n draws (Wilks one-sided tolerance interval)."""
    if k <= 0:
        return float(_beta.ppf(1 - delta, 1, n))
    if k >= n:
        return 1.0
    return float(_beta.ppf(1.0 - delta, k, n - k + 1))


def bound_threshold(controls: np.ndarray, fpr_nominal: float = 0.05,
                      delta: float = 0.05) -> dict:
    """Error-bound a data-chosen RNA-drop threshold at nominal FPR = fpr_nominal.

    threshold = k-th order statistic with k = round(fpr_nominal * n).
    Error-bounded ceiling is the Wilks order-statistic tolerance bound (sole lead).
    """
    s = np.sort(np.asarray(controls, float))
    n = len(s)
    k = int(round(fpr_nominal * n))
    k = max(1, min(k, n - 1))
    thr = float(s[k - 1])
    below = int((s <= thr).sum())
    return dict(n=n, k=k, threshold=thr, k_below=below, point_fpr=below / n,
                wilks_ceiling=wilks_tolerance_upper(k, n, delta))


def bound_threshold_at(controls: np.ndarray, fpr_ceiling: float = 0.05,
                           delta: float = 0.05) -> dict:
    """Least-conservative threshold whose Wilks-error-bounded FPR ceiling does not
    exceed `fpr_ceiling`.  Supersedes the preliminary point threshold: its FPR
    is bounded by `fpr_ceiling` with confidence 1 - delta, not merely estimated.
    """
    s = np.sort(np.asarray(controls, float))
    n = len(s)
    best_k = 0
    for k in range(1, n):
        if wilks_tolerance_upper(k, n, delta) <= fpr_ceiling:
            best_k = k
        else:
            break
    if best_k == 0:
        return dict(n=n, k=0, threshold=float(s[0]) - 1e-9, point_fpr=0.0,
                    wilks_ceiling=wilks_tolerance_upper(1, n, delta),
                    code_eligible=False)
    thr = float(s[best_k - 1])
    below = int((s <= thr).sum())
    return dict(n=n, k=best_k, threshold=thr, k_below=below, point_fpr=below / n,
                wilks_ceiling=wilks_tolerance_upper(best_k, n, delta),
                code_eligible=True)


# Correlation-aware widening (CP effective-sample)                            #

# Nearby variants are positionally correlated, so raw-count CP is optimistic.
# Primary widening: position-clustering (one independent draw per genomic
# position).  Stability check: block bootstrap over position + a design-effect
# (VIF) plug-in from the error-indicator autocorrelation.

def deff_from_acf(e_ordered: np.ndarray, max_lag: int = 5) -> float:
    """Design effect (variance-inflation factor) from the first `max_lag`
    autocorrelations:  DEFF = 1 + 2 * sum_l (1 - l/(L+1)) * rho_l  (Bartlett
    taper, truncated at the i.i.d. floor of 1)."""
    e = np.asarray(e_ordered, float)
    n = len(e)
    if n < 3 or e.std() == 0:
        return 1.0
    ec = e - e.mean()
    denom = (ec * ec).sum()
    if denom <= 0:
        return 1.0
    deff = 1.0
    L = min(max_lag, n - 1)
    for l in range(1, L + 1):
        rho = (ec[:-l] * ec[l:]).sum() / denom
        w = 1.0 - l / (L + 1.0)
        deff += 2.0 * w * rho
    return float(max(deff, 1.0))


def _cluster_cp(err: np.ndarray, pos: np.ndarray, delta: float) -> tuple:
    """Position-clustered CP: aggregate errors to a per-position mean, treat the
    number of distinct positions as the effective sample size.  Returns
    (n_clusters, cp_ceiling)."""
    per = pd.DataFrame({"pos": np.asarray(pos), "e": np.asarray(err, float)}
                       ).groupby("pos")["e"].mean()
    n_cl = len(per)
    k_cl = per.sum()
    return int(n_cl), ub_cp(k_cl, n_cl, delta)


def widen_group(err_ordered: np.ndarray, pos: np.ndarray, delta: float = 0.05,
                n_boot: int = 4000, block: int | None = None,
                seed: int = 20260710) -> dict:
    """Error-bounded CP ceiling under each correlation model (all CP-based).

    err_ordered : 0/1 errors ORDERED BY genomic position.
    pos         : matching genomic positions.
    Returns the independence-assuming reference plus the three widened ceilings,
    with `cp_ceiling` = the primary (position-clustered) value.
    """
    e = np.asarray(err_ordered, float)
    pos = np.asarray(pos)
    n = len(e)
    k = e.sum()
    out = dict(n=n, k=float(k), point=float(k / n) if n else np.nan)

    # independence-assuming reference (NOT the error-bounded ceiling)
    out["cp_raw_indep"] = ub_cp(k, n, delta)

    # (1) block bootstrap over position (circular blocks preserve local dep.)
    if n >= 10:
        rng = np.random.default_rng(seed)
        b = block or max(1, int(round(n ** (1.0 / 3.0))))
        nb = int(np.ceil(n / b))
        means = np.empty(n_boot)
        for i in range(n_boot):
            starts = rng.integers(0, n, size=nb)
            take = (starts[:, None] + np.arange(b)[None, :]).ravel() % n
            means[i] = e[take[:n]].mean()
        out["cp_block_boot"] = float(np.quantile(means, 1 - delta))
        out["block_len"] = b
    else:
        out["cp_block_boot"] = 1.0
        out["block_len"] = 0

    # (2) design-effect / VIF plug-in
    deff = deff_from_acf(e)
    n_eff = n / deff
    k_eff = out["point"] * n_eff
    out["deff"] = deff
    out["n_eff_vif"] = n_eff
    out["cp_vif"] = ub_cp(k_eff, n_eff, delta)

    # (3) position-clustering — primary
    n_cl, cp_cl = _cluster_cp(e, pos, delta)
    out["n_clusters"] = n_cl
    out["cp_cluster"] = cp_cl

    # the error-bounded ceiling is the primary (clustered) value
    out["cp_ceiling"] = cp_cl
    # stability spread across the three widened models
    widened = [out["cp_cluster"], out["cp_block_boot"], out["cp_vif"]]
    out["widen_spread"] = float(np.nanmax(widened) - np.nanmin(widened))
    return out


# Grid-error bounding engine with abstention                                   #

# Error-bounded ceiling for a call rate = position-clustered CP (correlation-aware).
# Abstain (return "uncertain") when n < MIN_N or the error-bounded ceiling is
# uninformative (> VACUOUS_CEILING).  Strata are ALWAYS within one gene and one
# variant class; genes and classes are never pooled.

VACUOUS_CEILING = 0.50
MIN_N = 10


def bound_group(err: np.ndarray, pos: np.ndarray, delta: float = 0.05,
                  vacuous: float = VACUOUS_CEILING, min_n: int = MIN_N) -> dict:
    """Error-bound one single-gene, single-class group given a 0/1 error vector
    `err` and matching genomic positions `pos`.  Error-bounded ceiling = CP on the
    position-clustered effective sample."""
    err = np.asarray(err, float)
    pos = np.asarray(pos)
    n = int(len(err))
    k = float(err.sum())
    if n == 0:
        return dict(n=0, k=0, point=np.nan, n_clusters=0, cp_ceiling=1.0,
                    abstain=True, abstain_reason="empty")
    point = k / n
    n_cl, cp_ceiling = _cluster_cp(err, pos, delta)
    abstain = (n < min_n) or (cp_ceiling > vacuous)
    reason = ("n<%d" % min_n) if n < min_n else (
        "vacuous(>%.2f)" % vacuous if cp_ceiling > vacuous else "")
    return dict(n=n, k=int(k), point=point, n_clusters=int(n_cl),
                cp_ceiling=float(cp_ceiling), abstain=bool(abstain),
                abstain_reason=reason)


def bound_grid(df: pd.DataFrame, group_cols, err_col: str, pos_col: str = "pos",
                 delta: float = 0.05, direction: str = "",
                 vacuous: float = VACUOUS_CEILING, min_n: int = MIN_N) -> pd.DataFrame:
    """Grid error-bound: split `df` by `group_cols` (which must include gene and a
    single variant class), bound the 0/1 `err_col` per cell using the
    position-clustered CP ceiling, and return one row per group.
    """
    rows = []
    for keys, sub in df.groupby(group_cols, dropna=False, observed=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        rec = dict(zip(group_cols, keys))
        rec["direction"] = direction
        rec.update(bound_group(sub[err_col].values, sub[pos_col].values,
                                  delta, vacuous, min_n))
        rows.append(rec)
    out = pd.DataFrame(rows)
    stat_cols = ["direction", "n", "k", "point", "n_clusters", "cp_ceiling",
                 "abstain", "abstain_reason"]
    return out[list(group_cols) + stat_cols].sort_values(
        list(group_cols)).reset_index(drop=True)
