# splice_downstream — finding splice-altering variants that DNA-based tools miss

A downstream, forward-looking application of the two-score method. It tests one
concrete claim, on rare variants in the splice region and exon core:

> Among variants that DNA-based practice would read as protein-coding
> (missense/synonymous), our screens' mRNA score flags a subset that lowers the
> **transcript**. For the ones an independent experimental database (SpliceVarDB)
> has data on, (1) how often is the RNA-drop flag right (**precision**), and
> (2) of the ones independently confirmed as splice-altering, how many does
> **SpliceAI miss** (**recovery**)?

The point: our two-score readout recovers real splice effects that sequence-based
practice drops — a validation of the method and a promising, additive use of it.

## Why this is the natural place to look

SpliceAI is weakest away from the canonical splice dinucleotides — inside the exon
core, where most of our candidates sit (1,383 of 1,614 are >3 nt from a junction).
SpliceVarDB is built for exactly this blind spot: most of its splice-altering
variants lie **outside** the canonical sites. So the overlap of "our RNA-drop
candidates" with "SpliceVarDB-confirmed splice-altering" is where a measurement-based
method should beat a sequence model.

## Pipeline

```
build_candidates.py     ../data/pooled_labeled.csv         -> data/candidate_variants.csv
   (offline)            missense/synonymous with rna_drop5 == True  (1,614 variants)

run_spliceai.py         data/candidate_variants.csv        -> data/spliceai_scores.csv
   (NEEDS NETWORK)      Broad SpliceAI-lookup API, hg38      ds_max per variant
                        default: only the SpliceVarDB-overlap set (~33 queries)

analyze.py              candidate + SpliceVarDB + SpliceAI -> results/*.csv, summary.json
   (offline)            precision + recovery
```

Run order:

```bash
python build_candidates.py      # rebuild candidates (offline)
python run_spliceai.py          # fetch SpliceAI for the overlap set (needs network)
python analyze.py               # precision + recovery (offline)
```

`run_spliceai.py` is the only step that needs the network. It is resumable
(skips variants already scored), retries transient 429/502 errors, and sleeps
between calls. Use `--all` to score all 1,614 candidates (only needed for the
full-candidate context number); the default scores just the SpliceVarDB-overlap
set, which is all the precision/recovery numbers require.

If the Broad API is down, any tool that emits a **`ds_max`** per variant works —
just produce `data/spliceai_scores.csv` with columns `key,gene,ds_max` (key =
`chr-pos-ref-alt`, hg38) and run `analyze.py`.

## Data

| file | what | source |
|---|---|---|
| `data/candidate_variants.csv` | 1,614 missense/synonymous variants with an mRNA drop | our pooled screens |
| `data/splicevardb_5genes.tsv` | SpliceVarDB export, 5 genes, hg38 | SpliceVarDB (splicevardb.org) |
| `data/spliceai_scores.csv` | SpliceAI ds_max per variant | Broad SpliceAI-lookup (you run) |

**Independence.** `analyze.py` drops any SpliceVarDB entry whose evidence DOI is one
of the five SGE screens this method is built on, so every confirmation is
independent of our own data. In the shipped export this removes 0 entries — none of
the screen DOIs appear (SpliceVarDB's evidence is classical minigene / RT-PCR /
MaPSy / MFASS / RNA-Seq assays) — but the filter is applied and reported regardless.

## SpliceAI thresholds

The published SpliceAI delta-score cutoffs (Jaganathan et al., Cell 2019) are 0.2
(high recall), 0.5 (recommended), 0.8 (high precision). "SpliceAI-negative" means
`ds_max` below the threshold. Recovery is reported at all three; the headline uses
0.2 — the most generous to SpliceAI, so a variant it misses at 0.2 it misses at any
clinical cutoff.

## Results (written by analyze.py)

| file | what |
|---|---|
| `results/overlap_by_gene.csv` | candidate × SpliceVarDB, gene × classification |
| `results/precision.csv` | precision of the RNA-drop flag (pooled + per gene) |
| `results/recovery.csv` | SpliceAI-negative fraction among confirmed, per threshold |
| `results/recovered_variants.csv` | the confirmed splice-altering **and** SpliceAI-negative variants |
| `results/summary.json` | headline numbers |
