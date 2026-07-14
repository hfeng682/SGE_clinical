# the reliability-map stage — AlphaGenome reliability map, tissue bridge, GTEx audit

Two pipelines, in order:

1. **Scoring** (`ag_common.py`, `pilot_run.py`, `run_shard.py`, `combine_shards.py`) —
   scores every variant in `../analysis/pooled_labeled.csv` (35,333) with
   AlphaGenome and produces one wide per-variant table (`alphagenome_per_variant.csv`):
   splice summaries plus tissue-resolved expression. **These four scripts are
   frozen; run them once to produce the table.**
2. **Analysis** (`analysis_common.py`, `reliability.py`, `tissue_bridge.py`,
   `gtex_audit.py`, `figures.py`) — everything downstream of that table: the
   reliability map, the reliability-gated tissue bridge, the GTEx
   double-counting audit, and the one figure. Outputs land in `results/`.

The narrative from start to end is `STEP4_main_story.md`. This README is the
run procedure.

## Install

```bash
pip install alphagenome pandas numpy
```

AlphaGenome is a hosted API (no GPU needed locally). Get a free non-commercial
key from Google DeepMind and keep it out of the scripts — pass it on the command
line or via `$ALPHAGENOME_API_KEY`. **The scripts never write the key anywhere.**

## What gets scored

Six scorers per variant, bundled into one request (frozen in `ag_common.py`):

| Scorer | Output | Purpose |
|---|---|---|
| `GeneMaskLFCScorer` | RNA_SEQ | expression log2 fold-change, tissue-resolved |
| `GeneMaskSplicingScorer` | SPLICE_SITES | splice-site class-probability change |
| `GeneMaskSplicingScorer` | SPLICE_SITE_USAGE | splice-site usage change |
| `SpliceJunctionScorer` | SPLICE_JUNCTIONS | junction-usage change |
| `GeneMaskActiveScorer` | RNA_SEQ | abundance (UTR/abundance branch) |
| `PolyadenylationScorer` | — | polyA (UTR branch) |

Each variant is scored in a 1 MB context window centred on its position
(covers every gene body here in full). AlphaGenome's tidy long output is reduced
to one wide row per variant by `ag_common.reduce_to_wide` — signed + magnitude
summaries per modality, a per-GTEx-tissue expression vector, and a
`has_gtex_track` flag for the later GTEx double-counting audit.

## Step 1 — pilot (~50 variants), do this first

```bash
python pilot_run.py --api-key "$ALPHAGENOME_API_KEY"
```

Writes `pilot_out/pilot_wide.csv` (the reduced table) and
`pilot_out/pilot_tidy.csv` (raw per-track scores). Send both back for a
sanity-check against the measured RNA score before the full run.

## Step 2 — full run, sharded across your API keys

You have ~5–6 keys. Run one process per key, each with a different `--shard-id`,
so each key scores ~1/N of the variants. Strided sharding keeps every gene
represented in every shard.

```bash
# 6 keys → 6 shards (~5,889 variants each). Run these in 6 terminals/processes:
python run_shard.py --api-key "$KEY0" --num-shards 6 --shard-id 0
python run_shard.py --api-key "$KEY1" --num-shards 6 --shard-id 1
python run_shard.py --api-key "$KEY2" --num-shards 6 --shard-id 2
python run_shard.py --api-key "$KEY3" --num-shards 6 --shard-id 3
python run_shard.py --api-key "$KEY4" --num-shards 6 --shard-id 4
python run_shard.py --api-key "$KEY5" --num-shards 6 --shard-id 5
```

- One progress line per variant; checkpoints every 50 variants.
- **Resumable:** rerun the same command and finished variants are skipped.
- Each shard writes `out/shard_XX_of_YY.csv`.
- `--seq-length {16KB,100KB,500KB,1MB}` trades speed for context (default 1MB).
- `--max-workers N` sets concurrent in-flight requests per key (default 4).

## Step 3 — combine into the end product

```bash
python combine_shards.py --num-shards 6
```

Merges all shards, checks every variant is covered exactly once, joins the
AlphaGenome columns onto `pooled_labeled.csv`, and writes
**`alphagenome_per_variant.csv`** — the single end-product wide table. Refuses to
write if variants are missing or errored (rerun those shards; pass
`--allow-partial` to override).

## the reliability-map stage — analysis (from the table to the end)

Once `alphagenome_per_variant.csv` exists, the analysis is four scripts run in
order (all read/write `results/`; no API, no network):

```bash
python reliability.py          # reliability map + missed demotions
python tissue_bridge.py  # reliability-gated tissue evidence for intronic variants
python gtex_audit.py     # GTEx provenance + collision audit
python figures.py        # the one figure (reads reliability map outputs)
```

- `reliability.py` runs a 2,000× bootstrap over ~40 strata (~4 min); the others are seconds.
- `analysis_common.py` holds every frozen definition (ground truth, predictors,
  distance bins, the reliability rule `AUROC CI-lower > 0.70`, disease tissues,
  the bundled GTEx v8 tissue panel). Change analysis choices there, in one place.
- `gtex_audit.py`'s variant-collision step reads a bundled eQTL cache
  (`results/gtex_eqtls.csv`, fetched from the GTEx v8 API during development);
  the provenance audit (the decisive finding) needs no external data.

## Files

**Scoring (frozen):**
- `ag_common.py` — shared: variant construction, scorer set, tidy→wide reduction
- `pilot_run.py` — ~50-variant pilot (not part of the full run)
- `run_shard.py` — sharded, resumable full-run scorer
- `combine_shards.py` — merge shards → `alphagenome_per_variant.csv`

**Analysis (from the table onward):**
- `analysis_common.py` — frozen definitions, metrics, reliability rule, tissue lists
- `reliability.py` — reliability map + missed demotions → `results/reliability_*.csv`, `reliability_summary.json`
- `tissue_bridge.py` — reliability-gated tissue evidence → `results/tissue_bridge.csv`
- `gtex_audit.py` — GTEx provenance + collision → `results/gtex_audit_summary.json`, `gtex_double_counting_audit.csv`
- `figures.py` — the one figure → `results/reliability_map.png`
- `STEP4_main_story.md` — the start-to-end narrative

## Coordinates (validated)

`pooled_labeled.csv` `pos` is 1-based and `ref`/`alt` are forward-strand — checked
34/34 against GRCh38 (Ensembl), including the four minus-strand genes. The only
transform is the `chrom → 'chr'` prefix, handled in `make_variant`.
