# validation/ — the four tests

Each test asks whether a guarantee is real. All reuse the frozen error-bound engine
(`model/`) unchanged and recompute their own tables. Run with `python run.py --validate`.

| script | the question | headline result |
|---|---|---|
| `bound_holds.py` (within-gene) | Does the 95% error bound actually cover the error on held-out controls? | **0.958** (23/24 strata) — holds. |
| `bound_holds.py` (leave-one-gene-out) | What happens if you break the design and bound a gene from *other* genes' controls? | **0.833** — under-covers on purpose. This is the **negative control** proving error rates are gene-specific, i.e. why the method never pools. |
| `tissue_discount.py` | Is the one-tier transfer discount the right size? | Tissue-invariant losses reproduce 87% vs tissue-variable 67–77%; the gap is **1.68 tiers (95% CI [−0.13, 3.74])**, which contains the one tier the method discounts. |
| `external_screens.py` | Does the survival axis transfer to screens the method never saw? | RAD51C 98.2% / 100%, DDX3X 100% / 79.4% agreement with independent ClinVar controls. |

**Reading the two coverage numbers together.** Within-gene (0.958) is how the method
is used and is the positive result. Leave-one-gene-out (0.833) is a *negative
control* — it is supposed to fail, and its failure is the evidence for the
no-pooling design. Neither is a checksum: both recompute error on controls the
thresholds never saw.

`external_data/` holds the two external one-score screens (RAD51C, DDX3X) and their
provenance.
