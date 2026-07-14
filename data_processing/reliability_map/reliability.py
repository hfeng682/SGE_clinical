#!/usr/bin/env python3
"""
reliability.py — AlphaGenome reliability map + missed-demotion count.

Ground truth: the measured mRNA-drop event (rna_drop5) on measurable exonic
variants. Predictor under test: AlphaGenome (its own splice scores; no second
model). We ask WHERE AG discriminates the measured event, resolved by variant
class and distance to the junction, and how many measured demotions AG would
miss at a fixed specificity.

Outputs (results/):
  reliability_discrimination.csv  per-(stratum x AG score) AUROC/AUPRC + bootstrap CI + reliable flag
  missed_demotions.csv      per-stratum missed-demotion counts at 90% / 95% specificity
  reliability_summary.json        frozen definitions, thresholds, overall AUROCs, reliable strata

Run:  python reliability.py
"""
import os
import sys
import json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analysis_common as C


def strata(m):
    """Yield (stratum_type, stratum_value, boolean mask) over the measurable frame."""
    yield ("overall", "all", np.ones(len(m), dtype=bool))
    for v in sorted(m["route_class"].dropna().unique()):
        yield ("route_class", v, (m["route_class"] == v).values)
    for v in sorted(m[C.DIST_BIN].dropna().unique()):
        yield ("dist_bin", v, (m[C.DIST_BIN] == v).values)
    # fine bins inside the exon core
    fb = C.fine_bin(m[C.DIST].values)
    for lab in C.FINE_LABELS:
        yield ("fine_bin", lab, (fb == lab))
    for v in sorted(m["gene"].dropna().unique()):
        yield ("gene", v, (m["gene"] == v).values)


def discrimination_table(m):
    m = m.copy()
    m["ag_splice_mag"] = C.ag_splice_mag(m)
    scores = {
        "ag_splice_mag": m["ag_splice_mag"].values,
        "splice_sites_maxabs": m["splice_sites_maxabs"].values,
        "splice_usage_maxabs": m["splice_usage_maxabs"].values,
        "splice_junctions_maxabs": m["splice_junctions_maxabs"].values,
        "expr_lfc_maxabs": m[C.EXPR_COL].values,
    }
    y = m[C.EVENT].values.astype(bool)
    rna_score = m["rna_score"].values if "rna_score" in m else np.full(len(m), np.nan)

    rows = []
    for stype, sval, mask in strata(m):
        n = int(mask.sum())
        ne = int(y[mask].sum())
        powered = (n >= C.MIN_N) and (ne >= C.MIN_EVENTS)
        # per-stratum reliability = best splice score's AUROC CI-lower > 0.70
        best_ci_lo = -np.inf
        for score_name, s in scores.items():
            if powered:
                a = C.auroc(y[mask], s[mask])
                lo, hi = C.bootstrap_ci(y[mask], s[mask])
                p = C.auprc(y[mask], s[mask])
                # Spearman of score vs continuous rna_score (sign check only)
                sm = ~np.isnan(rna_score[mask]) & ~np.isnan(s[mask])
                if sm.sum() > 10:
                    rr = pd.Series(s[mask][sm]).rank().values
                    rg = pd.Series(rna_score[mask][sm]).rank().values
                    rho = float(np.corrcoef(rr, rg)[0, 1])
                else:
                    rho = np.nan
            else:
                a = lo = hi = p = rho = np.nan
            if score_name in C.SPLICE_MAG_COLS + ["ag_splice_mag"] and not np.isnan(lo):
                best_ci_lo = max(best_ci_lo, lo)
            rows.append(dict(stratum_type=stype, stratum_value=sval, ag_score=score_name,
                             event=C.EVENT, n=n, n_events=ne,
                             event_rate=round(ne / n, 4) if n else np.nan,
                             auroc=a, auroc_lo=lo, auroc_hi=hi, auprc=p,
                             spearman_rna_score=rho, powered=powered))
        reliable = powered and (best_ci_lo > C.RELIABLE_CI_LO)
        for r in rows[-len(scores):]:
            r["reliable_flag"] = bool(reliable)
    return pd.DataFrame(rows)


def missed_demotions(m, disc):
    m = m.copy()
    m["ag_splice_mag"] = C.ag_splice_mag(m)
    y = m[C.EVENT].values.astype(bool)
    neg = m.loc[~y, "ag_splice_mag"].values
    thr = {spec: C.threshold_at_specificity(neg, spec) for spec in C.SPEC_TARGETS}

    rows = []
    for stype, sval, mask in strata(m):
        sub = m[mask]
        ys = sub[C.EVENT].values.astype(bool)
        dem = sub.loc[ys, "ag_splice_mag"].values  # demotions = true events
        nd = len(dem)
        if nd == 0:
            continue
        row = dict(stratum_type=stype, stratum_value=sval, n_demotions=nd)
        for spec in C.SPEC_TARGETS:
            missed = int((dem < thr[spec]).sum())
            row[f"missed_{int(spec*100)}spec"] = missed
            row[f"missed_rate_{int(spec*100)}spec"] = round(missed / nd, 4)
        rows.append(row)
    return pd.DataFrame(rows), thr


def main():
    os.makedirs(C.RESULTS_DIR, exist_ok=True)
    df = C.load_table()
    m = C.measurable(df)
    print(f"measurable exonic variants: {len(m)}  |  {C.EVENT} events: {int(m[C.EVENT].sum())}")

    disc = discrimination_table(m)
    disc.to_csv(os.path.join(C.RESULTS_DIR, "reliability_discrimination.csv"), index=False)

    miss, thr = missed_demotions(m, disc)
    miss.to_csv(os.path.join(C.RESULTS_DIR, "missed_demotions.csv"), index=False)

    # headline summary
    ov = disc[(disc.stratum_type == "overall") & (disc.event == C.EVENT)].set_index("ag_score")
    reliable_strata = sorted(
        f"{r.stratum_type}:{r.stratum_value}"
        for r in disc[disc.reliable_flag & (disc.ag_score == "ag_splice_mag")].itertuples())
    m2 = m.copy(); m2["ag_splice_mag"] = C.ag_splice_mag(m2)
    y = m2[C.EVENT].values.astype(bool)
    summary = dict(
        frozen_definitions=dict(
            ground_truth="rna_measured==True (measurable exonic); intronic variants are tissue-bridge targets, disjoint",
            n_measurable=int(len(m)), n_events=int(m[C.EVENT].sum()),
            event=C.EVENT, ag_splice_predictor="max(splice_sites,splice_usage,splice_junctions)_maxabs",
            ag_expression_predictor=C.EXPR_COL,
            auroc="Mann-Whitney rank statistic", auprc="trapezoidal PR",
            ci=f"percentile bootstrap {C.BOOTSTRAP_N}x seed {C.SEED}, 95%",
            power_filter=f"skip n<{C.MIN_N} or events<{C.MIN_EVENTS}",
            reliable_rule=f"best splice AUROC CI-lower > {C.RELIABLE_CI_LO}"),
        overall_auroc={k: float(ov.loc[k, "auroc"]) for k in ov.index},
        operating_thresholds={f"{int(s*100)}spec": thr[s] for s in C.SPEC_TARGETS},
        overall_sensitivity={f"{int(s*100)}spec": float((m2.loc[y, "ag_splice_mag"] >= thr[s]).mean())
                             for s in C.SPEC_TARGETS},
        reliable_strata=reliable_strata,
        n_strata_evaluated=int(disc[disc.powered & (disc.ag_score == "ag_splice_mag")].shape[0]))
    with open(os.path.join(C.RESULTS_DIR, "reliability_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # console report
    print(f"\noverall AUROC  splice={summary['overall_auroc']['ag_splice_mag']:.3f}  "
          f"expr={summary['overall_auroc']['expr_lfc_maxabs']:.3f}")
    g = disc[(disc.stratum_type == "dist_bin") & (disc.ag_score == "ag_splice_mag")]
    for r in g.itertuples():
        print(f"  dist_bin {r.stratum_value:26s} AUROC={r.auroc:.3f} "
              f"[{r.auroc_lo:.3f},{r.auroc_hi:.3f}]  reliable={r.reliable_flag}")
    fb = disc[(disc.stratum_type == "fine_bin") & (disc.ag_score == "ag_splice_mag")]
    print("  fine (exon core):", {r.stratum_value: round(r.auroc, 3) for r in fb.itertuples()})
    print(f"\nreliable strata (splice): {reliable_strata}")
    mid = miss[(miss.stratum_type == "fine_bin") & (miss.stratum_value == ">100")]
    if len(mid):
        r = mid.iloc[0]
        print(f"mid-exon >100nt demotions: {int(r.n_demotions)}, "
              f"missed@90spec={int(r.missed_90spec)} ({r.missed_rate_90spec:.1%})")


if __name__ == "__main__":
    main()
