# External screens — provenance

Two external one-score SGE screens used for the Axis 1 validation (the validation section). Neither is in the
five-gene development set; both were retrieved and processed in-session with the repo's own
frozen front-end.

## Sources (MaveDB, GRCh38)

| Gene | MaveDB score set | n variants | Transcript | Investigator call column(s) |
|---|---|---|---|---|
| RAD51C | `urn:mavedb:00000673-0-1` | 9,188 | ENST00000337432.9 | `functional_classification` |
| DDX3X | `urn:mavedb:00000658-0-1` ("SGE scores for clinical use") | 9,079 | ENST00000644876.2 | `SGE_prediction_of_variant_function_in_NDD_context` × `Confidence_of_functionally_abnormal_variant_prediction` |

RAD51C: Olvera-León et al., *Cell* 2024 (Sanger HAP1 SGE). DDX3X: X-linked, neurodevelopmental
disorder — a different disease area from the cancer-predisposition development set.

## How each file was produced

- `data/raw/<GENE>_scores.csv` — downloaded from the MaveDB API (`/score-sets/<urn>/scores`, CSV).
- `data/annotated/<GENE>_annotated.csv` — produced by the repo's VEP front-end
  (`../../data/annotation_pipeline/`, `transform()`), same code and Ensembl endpoint as the five
  development screens; a runtime `GENES` entry supplies transcript/chrom/strand (no source file
  edited). 9,188/9,188 and 9,079/9,079 annotated on the target transcript.
- `external_data/controls/external_clinical_controls.csv` — built by the frozen clinical-controls logic: ClinVar controls
  scanned from the repo's local `../../data/clinvar_gnomad/clinvar.vcf.gz` (whole-genome; the
  RAD51C/DDX3X windows are covered), gnomAD r4 allele frequencies via the variants connector,
  2★ P/LP-or-B/LB → `calibration`, everything else → `apply`.
- `data/label_maps.json` — the investigator-label → {Normal/Uncertain/LoF} mapping, with sources.
  The mapping direction was verified against each screen's score sign (RAD51C: negative =
  depleted = LoF; DDX3X: positive = abnormal).

## Scope

Only **Axis 1** (survival call → error-bounded ACMG strength) is exercised. These screens have **no
mmRNA-abundance score**, so **Axis 2** (mechanism route + tissue-transfer discount) is not
testable here — see `../README.md` the validation section. No genuinely external *two-score* screen exists yet
(`data/sge_two_score_landscape.csv`).

`validate_external_axis1.py` rebuilds the controls and recomputes all `results/` tables
from `data/raw/` + `data/annotated/` + the local ClinVar VCF.
