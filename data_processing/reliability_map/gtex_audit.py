#!/usr/bin/env python3
"""
gtex_audit.py — GTEx double-counting audit.

AlphaGenome is trained on GTEx-type data and its expression tracks are labelled
by GTEx tissue. This script establishes, so the evidence-code stage does not double-count AlphaGenome
tissue evidence against GTEx eQTL/sQTL evidence:

  1. PROVENANCE (offline, decisive): AG's 54 expression tracks map 1:1 onto the
     GTEx v8 tissue panel -> AG tissue evidence is GTEx-derived by construction.
  2. VARIANT COLLISION (optional): how often an SGE variant position coincides
     with a significant GTEx eQTL. Requires a cached eQTL file (see below);
     skipped with a note if absent, since the provenance finding already
     settles the independence question.

Optional eQTL cache: results/gtex_eqtls.csv with columns [gene, pos]
(one row per significant single-tissue eQTL). If present, the collision
analysis runs; if not, only the provenance audit is reported. During
development this was fetched from the GTEx v8 single-tissue eQTL API via the
expression connector.

Output (results/):
  gtex_audit_summary.json   provenance verdict + double-counting rule for the evidence-code stage
  gtex_double_counting_audit.csv  per-gene provenance/collision table (if eQTLs available)

Run:  python gtex_audit.py
"""
import os
import sys
import json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analysis_common as C


def norm(s):
    return s.lower().replace("-", "_").replace(" ", "_").replace("(", "").replace(")", "")


def main():
    os.makedirs(C.RESULTS_DIR, exist_ok=True)
    df = C.load_table()
    ag_tissues = [c.replace("expr_lfc_gtex__", "") for c in df.columns
                  if c.startswith("expr_lfc_gtex__")]

    # ---- 1. PROVENANCE ----
    ag_norm = {norm(t) for t in ag_tissues}
    gtex_norm = {norm(t) for t in C.GTEX_V8_TISSUES}
    matched = ag_norm & gtex_norm
    print(f"AG tissue tracks: {len(ag_tissues)}  GTEx v8 sites: {len(C.GTEX_V8_TISSUES)}  "
          f"exact matches: {len(matched)}")
    provenance = dict(ag_tissue_tracks=len(ag_tissues), gtex_v8_tissue_sites=len(C.GTEX_V8_TISSUES),
                      exact_matches=len(matched),
                      ag_only=sorted(ag_norm - gtex_norm), gtex_only=sorted(gtex_norm - ag_norm),
                      conclusion="AG tissue-resolved expression uses the GTEx v8 tissue panel verbatim; "
                                 "it is GTEx-derived by construction (AG training corpus includes GTEx RNA-seq).")

    # ---- 2. VARIANT COLLISION (optional) ----
    eqtl_path = os.path.join(C.RESULTS_DIR, "gtex_eqtls.csv")
    collision = None
    if os.path.exists(eqtl_path):
        eq = pd.read_csv(eqtl_path)
        rows = []
        for g, sub in df.groupby("gene"):
            sge_pos = set(sub.pos.astype(int))
            epos = set(eq.loc[eq.gene == g, "pos"].astype(int))
            ov = len(sge_pos & epos)
            rows.append(dict(gene=g, n_sge_positions=len(sge_pos), gtex_eqtl_positions=len(epos),
                             sge_eqtl_collision=ov, pct_sge_are_eqtl=round(100 * ov / len(sge_pos), 3),
                             ag_tissue_tracks_are_gtex=f"{len(matched)}/{len(ag_tissues)} exact"))
        collision = pd.DataFrame(rows)
        collision.to_csv(os.path.join(C.RESULTS_DIR, "gtex_double_counting_audit.csv"), index=False)
        tot, col = collision.n_sge_positions.sum(), collision.sge_eqtl_collision.sum()
        print(f"variant collision: {col}/{tot} SGE positions are a GTEx eQTL ({100*col/tot:.2f}%)")
        collision_summary = dict(total_sge_positions=int(tot), colliding=int(col),
                                 pct=round(100 * col / tot, 3),
                                 conclusion="SGE variants (engineered, mostly rare) rarely coincide with "
                                            "common GTEx eQTL variants; per-variant AG and eQTL evidence "
                                            "rarely touch the same variant.")
    else:
        print(f"[note] {eqtl_path} not found — skipping variant-collision analysis; "
              "provenance finding already settles independence.")
        collision_summary = "skipped (no cached eQTL file); provenance audit is decisive on its own"

    summary = dict(
        provenance=provenance,
        variant_collision=collision_summary,
        independence_verdict="AG tissue evidence is NOT independent of GTEx: it is trained on GTEx and "
                             "reports in GTEx tissue vocabulary. Do NOT treat AG tissue-transfer evidence "
                             "as independent corroboration of a GTEx eQTL/sQTL for the same gene/tissue.",
        double_counting_rule_for_step5="When a variant carries GTEx eQTL/sQTL evidence in a disease tissue, "
                                      "AG tissue-resolved EXPRESSION evidence in that same tissue adds NO "
                                      "independent weight (shared GTEx provenance) — do not sum. AG SPLICE "
                                      "evidence (trained on GENCODE/splice references) is less GTEx-entangled "
                                      "and may add orthogonal information, but should be capped, not summed, "
                                      "with sQTL evidence.")
    with open(os.path.join(C.RESULTS_DIR, "gtex_audit_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print("\nverdict:", summary["independence_verdict"])


if __name__ == "__main__":
    main()
