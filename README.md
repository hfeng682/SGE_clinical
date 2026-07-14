# Two-score SGE → mechanism-labelled clinical evidence

Turn two-score saturation genome editing (SGE) screens into ACMG/AMP
functional-evidence codes (PS3 / BS3) that carry a **mechanism label**
(protein vs RNA) and a **finite-sample error bound**.


## Please read the manuscript before running this repo.

```
python run.py              # error bounds → strength → evidence codes  (stages 2–4)
python run.py --all        # also rebuild the clinical controls from ClinVar/gnomAD (stage 1)
python run.py --validate   # also run the four validation tests
```

The single output of record is
`results/evidence_codes/evidence_codes.csv` — one code per held-out variant.
Every table and figure in the manuscript regenerates from `results/`.

## Repository layout (organized by function)

```
run.py                  single entry point; chains the stages below

data/                   the inputs, frozen once
  pooled_labeled.csv        five screens, PASS-filtered, VEP-annotated, one shared call
  clinical_controls.csv     ClinVar 2★ / gnomAD controls + the frozen calibration/apply split
  alphagenome_scores.csv    per-variant AlphaGenome splice/expression predictions
  clinvar_gnomad/           raw ClinVar VCF + gnomAD r4 (to rebuild the controls)
  frozen_characterization/  derivation records not on the critical path

data_processing/        build the inputs from raw data
  build_clinical_controls.py   join screens to ClinVar/gnomAD, apply the frozen split
  summarize_controls.py        QC summary of the control set
  annotation/                  raw screen → VEP consequence → distance-to-junction
  reliability_map/             AlphaGenome reliability map (where the predictor is trusted)

model/                  the method
  bounds.py                finite-sample error bound on the survival call (Clopper–Pearson)
  strength.py              bounded error rate → ACMG strength via OddsPath
  mechanism.py             mechanism route (protein vs RNA) from the two scores
  tissue.py                tissue-transfer rule (one-tier discount for RNA route)
  build_bounds.py          → results/bounds/
  build_strength.py        → results/strength/
  build_evidence_codes.py  → results/evidence_codes/   (the final codes)

validation/             the four tests
  bound_holds.py           does the error bound hold on held-out controls?
                           (within-gene + leave-one-gene-out negative control)
  tissue_discount.py       is the one-tier transfer discount the right size?
  external_screens.py      does the survival axis transfer to screens never seen?
  external_data/           two external one-score screens (RAD51C, DDX3X)

utils/                  shared helpers
  paths.py                 single source of truth for every file path
  variant_keys.py          the variant key (gene, chrom, pos, ref, alt)

results/                regenerated outputs (bounds, strength, evidence_codes, reliability_map, validation)
manuscript/             MANUSCRIPT.md, figures/ (2), tables/ (11)
docs/                   methods notes (de-circularization notes, verification report)
```

