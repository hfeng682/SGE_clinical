"""
build_bounds.py -- driver: the finite-sample error bounds on the survival call.

Uses the intrinsic controls already in the pooled screen:
  * synonymous variants are functional negatives  -> a false-LoF call is an error
  * nonsense variants are functional positives    -> a false-Normal call is an error
For each gene x class (and, at the finer level, x distance bin) it bounds that
binary call-error rate with the position-clustered Clopper-Pearson ceiling
(bounds.bound_group / bound_grid), and reports the correlation-model
stability of the ceiling (bounds.widen_group).

Outputs (results/bounds/):
  error_bounds_by_distance.csv   gene x class x dist_bin x direction (fine grid)
  error_bounds_per_gene.csv      gene x class x direction (pooled over distance)
  correlation_stability.csv      the ceiling under 4 correlation models

NB the RNA-drop threshold and its Wilks tolerance bound were frozen upstream
(the cut lives in tissue.RNA_THR5 and the flag `rna_drop5` in the pooled data);
its derivation table is carried in data/frozen_characterization/rna_drop_threshold.csv.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))          # model/ (bounds.py)
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))       # repo root (utils)
from utils import paths
import bounds as B

# the two intrinsic control classes and the error each defines
NEG_CLASS, NEG_DIR = "synonymous", "false_LoF"     # negatives: a LoF call is wrong
POS_CLASS, POS_DIR = "nonsense",  "false_Normal"   # positives: a Normal call is wrong


def _err(sub: pd.DataFrame, wrong_call: str) -> np.ndarray:
    return (sub["call"] == wrong_call).astype(float).values


def per_gene_rollup(pool: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene, cls, direction, wrong in [
        (None, NEG_CLASS, NEG_DIR, "LoF"), (None, POS_CLASS, POS_DIR, "Normal")]:
        for g, sub in pool[pool.route_class == cls].groupby("gene"):
            rec = dict(level="L1_rollup_per_gene", gene=g, route_class=cls, direction=direction)
            rec.update(B.bound_group(_err(sub, wrong), sub.pos.values))
            rows.append(rec)
    cols = ["level","gene","route_class","direction","n","k","point","n_clusters",
            "cp_ceiling","abstain","abstain_reason"]
    return pd.DataFrame(rows)[cols]


def by_distance_grid(pool: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cls, direction, wrong in [(NEG_CLASS, NEG_DIR, "LoF"), (POS_CLASS, POS_DIR, "Normal")]:
        sub = pool[pool.route_class == cls].copy()
        sub["err"] = _err(sub, wrong)
        grid = B.bound_grid(sub, ["gene","route_class","dist_bin"], "err", "pos",
                              direction=direction)
        grid.insert(0, "level", "L0_grid")
        rows.append(grid)
    out = pd.concat(rows, ignore_index=True)
    return out


def correlation_stability(pool: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cls, direction, wrong in [(NEG_CLASS, NEG_DIR, "LoF"), (POS_CLASS, POS_DIR, "Normal")]:
        for g, sub in pool[pool.route_class == cls].groupby("gene"):
            sub = sub.sort_values("pos")
            w = B.widen_group(_err(sub, wrong), sub.pos.values)
            rows.append(dict(gene=g, route_class=cls, direction=direction, **w))
    df = pd.DataFrame(rows)
    keep = ["gene","route_class","direction","n","k","point","cp_raw_indep",
            "cp_block_boot","deff","cp_vif","n_clusters","cp_cluster","cp_ceiling","widen_spread"]
    return df[keep]


def main() -> None:
    pool = pd.read_csv(paths.POOLED)
    paths.RESULTS_BOUNDS.mkdir(parents=True, exist_ok=True)
    per_gene_rollup(pool).to_csv(paths.RESULTS_BOUNDS / "error_bounds_per_gene.csv", index=False)
    by_distance_grid(pool).to_csv(paths.RESULTS_BOUNDS / "error_bounds_by_distance.csv", index=False)
    correlation_stability(pool).to_csv(paths.RESULTS_BOUNDS / "correlation_stability.csv", index=False)
    print("wrote error_bounds_per_gene.csv, error_bounds_by_distance.csv, correlation_stability.csv")


if __name__ == "__main__":
    main()
