# Two-score SGE → mechanism-labelled clinical evidence

Turn two-score saturation genome editing (SGE) screens into ACMG/AMP
functional-evidence codes (PS3 / BS3) that carry a **mechanism label**
(protein vs RNA) and a **finite-sample error bound**.


## What it does, in one paragraph

Each variant in an SGE screen comes with two scores: a survival/function score
and an mRNA-abundance score. From the survival score the method makes a
three-class call (normal / uncertain / loss-of-function) and puts a
distribution-free error bound on how often that call is wrong, using clinical
controls of the variant's own class, then maps the bounded rate onto the ACMG
strength ladder via OddsPath. From the two scores together it reads a
**mechanism route** — a loss acting through a broken protein (transfers to any
tissue) or through disrupted RNA (tissue-variable) — and down-weights an
RNA-route call by one ACMG tier unless it is corroborated in disease-relevant
tissue. The output is one mechanism-labelled, transfer-adjusted evidence code
per variant.

## Run it

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

## Design principles

- **Bound per gene, never pool.** Error rates are gene-specific; the
  leave-one-gene-out test (`validation/bound_holds.py`) is the negative control
  that shows pooling breaks the guarantee.
- **The call uses the assay only.** The three-class call comes from the survival
  score and the investigators' published thresholds; clinical labels enter only
  downstream (control selection, strength calibration), so no variant grades
  itself.
- **Abstain rather than over-claim.** A class earns a code only if its bounded
  error rate clears an OddsPath rung; otherwise the method returns *uncertain*.
- **One command, frozen inputs.** `run.py` regenerates every result of record
  byte-for-byte — a software-integrity guarantee, kept separate from the
  empirical validations in `validation/`.
