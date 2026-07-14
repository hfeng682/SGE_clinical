# data_processing/ — build the inputs from raw data

Turns raw screens and public databases into the frozen files in `data/`. You only
need this to rebuild inputs from scratch (`python run.py --all`); the shipped
`data/` files are already built.

| script / folder | what it does | writes |
|---|---|---|
| `annotation/` | raw screen → normalized coordinates → VEP molecular consequence → distance-to-junction bin | `data/annotated_screens/` |
| `build_clinical_controls.py` | join assayed variants to ClinVar 2★ / gnomAD r4, assign control roles, apply the frozen calibration/apply split | `data/clinical_controls.csv` |
| `summarize_controls.py` | QC summary of the control set (counts, assay-vs-clinical agreement) | `results/control_summary.png` |
| `reliability_map/` | build the AlphaGenome **reliability map** — where, by distance to junction, the splice predictor discriminates well enough (lower AUROC bound > 0.70) to be trusted for tissue corroboration | `results/reliability_map/` |

**Reliability map, in plain terms.** AlphaGenome is only allowed to corroborate an
RNA-route loss where it has been shown to work — very close to the splice junction.
The map measures its discrimination against the screens' own mRNA readout, by
distance, and gates corroboration accordingly. The scoring-side scripts
(`reliability_map/ag_common.py`, `run_shard.py`, …) need the AlphaGenome API and are
kept for provenance; the analysis reads the cached scores in `data/`.
