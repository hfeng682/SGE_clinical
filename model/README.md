# model/ — the method

The four steps that turn one variant's two scores into one evidence code. Each
`build_*.py` runs a stage and writes to `results/`; the plain modules hold the logic.

| module | the step | what it computes |
|---|---|---|
| `bounds.py` | **error bound** | Distribution-free, finite-sample upper bound (exact-binomial Clopper–Pearson) on how often the survival call is wrong, on a **position-clustered** effective sample so correlated nearby variants don't inflate the count. |
| `strength.py` | **strength** | Maps the bounded error rate onto the ACMG strength ladder (Supporting / Moderate / Strong) via ClinGen's OddsPath thresholds. Abstains if no rung is cleared. |
| `mechanism.py` | **mechanism route** | Reads the two scores + molecular consequence into a route: protein (missense, in-frame; or nonsense → NMD, which is constitutive) vs RNA (splice/abundance disruption). |
| `tissue.py` | **tissue transfer** | The transfer rule: a protein-route loss keeps its tier; an RNA-route loss is down-weighted **one ACMG tier** unless corroborated in disease-relevant tissue (a second measurement condition, or a reliable AlphaGenome call). |

**Build scripts**

| script | writes |
|---|---|
| `build_bounds.py` | `results/bounds/` — per-gene and by-distance error-rate ceilings |
| `build_strength.py` | `results/strength/` — the per-class OddsPath codes |
| `build_evidence_codes.py` | `results/evidence_codes/evidence_codes.csv` — **the final output**: merge the two scores, run mechanism → strength → tissue, emit one code per variant |

**What is deliberately *not* here.** The RNA-route conjunction (both scores must
agree) is used only to route a loss and set the tissue-transfer rule — never to
tighten a strength tier. Requiring both scores to agree can only remove false
RNA-route positives, so the published tiers rest on the survival axis alone.
