#!/usr/bin/env python3
"""
run.py -- one entry point for the whole method.

From the frozen inputs in data/, regenerate every result of record:

    stage                     module                          reads                         writes
    ----------------------------------------------------------------------------------------------------------
    1  clinical controls      data_processing/build_clinical_controls.py   data/clinvar_gnomad, screens   data/clinical_controls.csv
    2  error bounds           model/build_bounds.py           data/pooled_labeled.csv        results/bounds/
    3  OddsPath strength      model/build_strength.py         data/clinical_controls.csv     results/strength/
    4  evidence codes         model/build_evidence_codes.py   pooled + controls + reliability results/evidence_codes/

The AlphaGenome reliability map (results/reliability_map/) is produced separately
(model scoring needs the AlphaGenome API); its frozen outputs ship in the repo and
stage 4 reads them. See data_processing/reliability_map/RUN_PROCEDURE.md.

Usage
    python run.py              # run stages 2-4 (controls already built & frozen)
    python run.py --all        # also rebuild clinical controls from ClinVar/gnomAD (stage 1)
    python run.py --validate   # run stages 2-4 then the out-of-sample checks in validation/
"""
from __future__ import annotations
import sys, subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent


def _run(label: str, script: str) -> None:
    print(f"\n=== {label} ===  ({script})")
    r = subprocess.run([sys.executable, str(REPO / script)], cwd=str(REPO))
    if r.returncode != 0:
        raise SystemExit(f"stage failed: {script} (exit {r.returncode})")


def main(argv) -> None:
    do_all = "--all" in argv
    do_val = "--validate" in argv

    if do_all:
        _run("STAGE 1  clinical controls", "data_processing/build_clinical_controls.py")

    _run("STAGE 2  error bounds",        "model/build_bounds.py")
    _run("STAGE 3  OddsPath strength",   "model/build_strength.py")
    _run("STAGE 4  evidence codes",      "model/build_evidence_codes.py")

    if do_val:
        # bound_holds.py runs BOTH the within-gene bound-holds check and the
        # leave-one-gene-out (across-genes) check.
        for lbl, s in [("does the bound hold (within-gene + across-genes)", "validation/bound_holds.py"),
                       ("tissue discount", "validation/tissue_discount.py"),
                       ("external screens", "validation/external_screens.py")]:
            _run(f"VALIDATION  {lbl}", s)

    print("\ndone. final table: results/evidence_codes/evidence_codes.csv")


if __name__ == "__main__":
    main(sys.argv[1:])
