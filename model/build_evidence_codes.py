"""
build_evidence_codes.py -- assemble the final per-variant evidence codes.

Runs the five stages over the held-out reclassification variants and writes the
deliverable table. This module assembles already-bounded quantities; it
introduces no new statistics. Every number it emits is traceable to a frozen
upstream source:

  * the three-class survival call and the mechanism flags -> data/pooled_labeled.csv
  * the clinical role (calibration vs apply) and 2* labels -> data/clinical_controls.csv
  * the class-matched error-bounded OddsPath (build_strength) -> results/strength/oddspath_calls_per_class.csv
  * the whole-gene OddsPath fallback (build_strength)         -> results/strength/oddspath_per_gene.csv
  * the reliability-gated AG tissue evidence                  -> results/reliability_map/tissue_bridge.csv

Stages
  STAGE 1  mechanism       : each variant -> a mechanism           (mechanism.route_mechanism)
  STAGE 2  strength ladder : mechanism + gene x class -> ACMG tier (mechanism.assign_strength)
  STAGE 3  SVI code        : mechanism + direction -> SVI code     (mechanism.svi_code)
  STAGE 4  tissue discount : RNA-route losses down-weighted unless corroborated (tissue.py)
  STAGE 5  abstain         : anything the method cannot resolve -> uncertain


Design decisions
  * Dual-mechanism (protein-visible LoF that also drops mRNA -- the 501 cohort):
    the functional tier is kept (both scores agree on loss), and only the
    transfer tier is subject to the discount. Mechanism does not erase the
    protein evidence; it re-grades transferability.
  * Tissue discount magnitude: exactly one ACMG tier for an uncorroborated
    RNA-route loss. A protein-route loss is never discounted.
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))          # model/
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))       # repo root (utils)
from utils import paths
from utils.variant_keys import KEY, key as _key
import mechanism as M
import tissue as T


def load_merged() -> pd.DataFrame:
    """Merge the frozen per-variant sources on the unique variant key."""
    pool = pd.read_csv(paths.POOLED)
    ctrl = pd.read_csv(paths.CLINICAL_CONTROLS,
                       usecols=KEY + ["clin_label", "gold_stars", "role", "calib_fold"])
    pool["key"], ctrl["key"] = _key(pool), _key(ctrl)
    df = pool.merge(ctrl[["key", "clin_label", "gold_stars", "role", "calib_fold"]],
                    on="key", how="left")
    tb = pd.read_csv(paths.RESULTS_RELIABILITY / "tissue_bridge.csv")
    tb["key"] = _key(tb)
    df = df.merge(tb[["key", "ag_reliable", "ag_splice_disruptor", "transfer_verdict",
                      "disease_expr_lfc_mean", "tissue_lfc_std"]], on="key", how="left")
    assert df.role.notna().all(), "every variant must carry a clinical role"
    return df


def assemble(df: pd.DataFrame) -> pd.DataFrame:
    """Run Stages 1-5 and return the per-variant evidence-code frame."""
    per_class, per_gene = M.load_strength_sources()

    # Stage 1: mechanism
    df = df.copy()
    df["route"] = df.apply(M.route_mechanism, axis=1)
    df["is_lof_call"] = df["route"].map(M.ROUTE_IS_LOF)
    df["is_rna_route"] = df["route"].map(M.ROUTE_IS_RNA)
    df["direction"] = np.where(df["route"] == "benign", "Normal",
                               np.where(df["is_lof_call"], "LoF", "none"))

    recs = []
    for _rt in df.itertuples(index=False):
        r = _rt._asdict()
        route, direction = r["route"], r["direction"]
        rec = dict(functional_tier="uncertain", oddspath_lo=np.nan,
                   strength_source="n/a", svi_code="none",
                   tissue_corroborated=False, corroboration_source="none",
                   transfer_tier="uncertain", tier_changed_by_tissue=False)

        if route == "uncertain":
            rec["strength_source"] = "guaranteed_uncertain"
            recs.append(rec); continue

        # Stage 2: error-bounded strength (functional tier)
        tier, oplo, src = M.assign_strength(r["gene"], r["route_class"], direction, route,
                                            per_class, per_gene)
        rec["functional_tier"], rec["oddspath_lo"], rec["strength_source"] = tier, oplo, src
        # Stage 3: SVI code
        rec["svi_code"] = M.svi_code(route, direction)

        # Stage 4: tissue-transfer discount
        if route == "benign":
            rec["transfer_tier"] = tier                 # a normal call transfers everywhere
        elif not r["is_rna_route"]:
            rec["transfer_tier"] = tier                 # protein route: tissue-invariant
        else:
            corr, csrc = T.tissue_corroboration(r)       # RNA route (incl. dual)
            rec["tissue_corroborated"], rec["corroboration_source"] = corr, csrc
            if corr:
                rec["transfer_tier"] = tier
            else:
                rec["transfer_tier"] = M.downweight(tier, 1)
                rec["tier_changed_by_tissue"] = (rec["transfer_tier"] != tier)
        recs.append(rec)

    add = pd.DataFrame(recs, index=df.index)
    out = pd.concat([df, add], axis=1)

    def final_code(row):
        if row["route"] == "uncertain" or row["transfer_tier"] == "uncertain":
            return "uncertain"
        return f"{row['svi_code']}_{row['transfer_tier']}"
    out["final_code"] = out.apply(final_code, axis=1)
    return out


# MECHANISM-BLIND COMPARISON                                                   #
# The mechanism-blind tier is this same method with the second (mRNA) score    #
# switched off: it keeps the survival call + OddsPath tier but never sees the  #
# RNA score, so it cannot detect the dual cohort and never applies the tissue  #
# discount. The count of variants whose tier differs is therefore, by          #
# definition, the number of RNA-route losses the tissue rule down-weighted --  #
# a decomposition of the discount, not an independent finding.                 #
def mechanism_blind_tier(row) -> str:
    if row["route"] == "uncertain" or not (row["is_lof_call"] or row["route"] == "benign"):
        return "uncertain"
    return row["functional_tier"]


def build_mechanism_vs_blind(out: pd.DataFrame) -> pd.DataFrame:
    o = out.copy()
    o["blind_tier"] = o.apply(mechanism_blind_tier, axis=1)
    o["aware_tier"] = o["transfer_tier"]
    o["tier_delta"] = (o["aware_tier"].map(M.TIER_RANK) - o["blind_tier"].map(M.TIER_RANK))
    o["mechanism_changed_strength"] = o["tier_delta"] != 0
    return o


def emit_summary(apply: pd.DataFrame, mvb: pd.DataFrame) -> dict:
    c501 = apply[(apply.route == "dual") & (apply.route_class == "missense")]
    return {
        "n_apply": int(len(apply)),
        "route_counts": apply.route.value_counts().to_dict(),
        "final_code_counts": apply.final_code.value_counts().to_dict(),
        "strength_source_counts": apply.strength_source.value_counts().to_dict(),
        "corroboration_source_counts":
            apply[apply.is_rna_route].corroboration_source.value_counts().to_dict(),
        "n_mechanism_changed": int(mvb.mechanism_changed_strength.sum()),
        "mechanism_changed_by_route":
            mvb[mvb.mechanism_changed_strength].route.value_counts().to_dict(),
        "mechanism_changed_by_gene":
            mvb[mvb.mechanism_changed_strength].gene.value_counts().to_dict(),
        "cohort501_apply": int(len(c501)),
        "cohort501_downgraded": int((c501.transfer_tier.map(M.TIER_RANK)
                                     < c501.functional_tier.map(M.TIER_RANK)).sum()),
        "cohort501_corroborated": int(c501.tissue_corroborated.sum()),
    }


def main() -> None:
    df = load_merged()
    out = assemble(df)
    apply = out[out.role == "apply"].copy()          # frozen held-out discipline
    assert (apply.role == "apply").all(), "apply table leaked a non-apply row"
    assert apply.calib_fold.isna().all(), "apply table leaked a calibration-fold variant"

    prov_cols = (KEY + ["gene", "route_class", "consequence", "call",
                 "protein_altering", "rna_measured", "rna_drop5",
                 "route", "direction", "is_lof_call", "is_rna_route",
                 "strength_source", "oddspath_lo", "functional_tier",
                 "svi_code",
                 "tissue_corroborated", "corroboration_source",
                 "transfer_tier", "tier_changed_by_tissue", "final_code"])
    prov_cols = [c for c in prov_cols if c in apply.columns]
    seen, ordered = set(), []
    for c in prov_cols:
        if c not in seen:
            ordered.append(c); seen.add(c)
    paths.RESULTS_CODES.mkdir(parents=True, exist_ok=True)
    apply[ordered].to_csv(paths.RESULTS_CODES / "evidence_codes.csv", index=False)

    mvb = build_mechanism_vs_blind(apply)
    mvb_cols = ordered + ["blind_tier", "aware_tier", "tier_delta",
                          "mechanism_changed_strength"]
    mvb_cols = [c for c in mvb_cols if c in mvb.columns]
    mvb[mvb_cols].to_csv(paths.RESULTS_CODES / "mechanism_vs_blind.csv", index=False)

    summary = emit_summary(apply, mvb)
    with open(paths.RESULTS_CODES / "evidence_codes_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print("evidence_codes.csv :", (apply.role == "apply").sum(), "apply variants")
    print("mechanism dist:", apply.route.value_counts().to_dict())
    print("final_code dist:", apply.final_code.value_counts().to_dict())
    print("mechanism changed strength:", int(mvb.mechanism_changed_strength.sum()))


if __name__ == "__main__":
    main()
