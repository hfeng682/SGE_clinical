"""
ag_common.py — shared helpers for AlphaGenome scoring of the pooled SGE variants.

Used by pilot_run.py, run_shard.py, and combine_shards.py so the variant
construction, the frozen scorer set, and the tidy->wide reduction ("how we
summarise AlphaGenome's output per variant") all live in exactly one place.

Nothing here connects to the network or holds an API key; the client is created
in the run scripts from the key the user passes on the command line.

Reduction design (per variant, keyed by `accession`):
  Expression (GeneMaskLFCScorer, RNA_SEQ, signed log2 fold-change ALT vs REF):
    - summary scalars across all RNA-seq tracks for the target gene
    - a tissue-resolved vector: one column per GTEx tissue (mean over its tracks)
  Expression "active" (GeneMaskActiveScorer, RNA_SEQ): magnitude summaries
  Splicing (GeneMaskSplicingScorer): SPLICE_SITES and SPLICE_SITE_USAGE magnitudes
  Splice junctions (SpliceJunctionScorer): junction-usage change magnitude
  Polyadenylation (PolyadenylationScorer): summaries for the UTR/abundance branch
  Provenance: count/flag of GTEx-derived tracks (for the later GTEx double-counting audit)

Coordinate note (validated 34/34 against GRCh38 via Ensembl): pooled_labeled.csv
`pos` is 1-based and `ref`/`alt` are on the forward genomic strand, which is
exactly what genome.Variant expects. The only transform is chrom -> 'chr' prefix.
"""
from __future__ import annotations

import time
import numpy as np
import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

# Configuration

# Sequence context window centred on the variant. Genes here span <=81 kb, so
# 1 MB (the model's full native context) comfortably covers the gene body plus
# regulatory flanks for every variant. Downshift with --seq-length for speed:
#   16KB / 100KB / 500KB / 1MB.
SEQUENCE_LENGTHS = {
    "16KB": dna_client.SEQUENCE_LENGTH_16KB,
    "100KB": dna_client.SEQUENCE_LENGTH_100KB,
    "500KB": dna_client.SEQUENCE_LENGTH_500KB,
    "1MB": dna_client.SEQUENCE_LENGTH_1MB,
}
DEFAULT_SEQ_LENGTH = "1MB"

# The frozen project scorer set. We deliberately do NOT run the full 19
# "recommended" scorers (ATAC/DNASE/ChIP/CAGE/contact maps are irrelevant to a
# splice/expression mechanism question). All six fit in one request
# (MAX_VARIANT_SCORERS_PER_REQUEST = 20).
SCORER_KEYS = [
    "RNA_SEQ",           # GeneMaskLFCScorer     -> expression log2FC, tissue-resolved
    "SPLICE_SITES",      # GeneMaskSplicingScorer -> splice-site class prob change
    "SPLICE_SITE_USAGE", # GeneMaskSplicingScorer -> splice-site usage change
    "SPLICE_JUNCTIONS",  # SpliceJunctionScorer   -> junction usage change
    "RNA_SEQ_ACTIVE",    # GeneMaskActiveScorer   -> abundance (UTR/abundance branch)
    "POLYADENYLATION",   # PolyadenylationScorer  -> polyA (UTR branch)
]

# Columns we read out of pooled_labeled.csv to build variants.
INPUT_KEY_COLS = ["accession", "chrom", "pos", "ref", "alt", "gene"]


def get_scorers():
    """Return the frozen list of AlphaGenome variant scorers for this project."""
    rec = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    return [rec[k] for k in SCORER_KEYS]


def load_variants(input_csv):
    """Load pooled_labeled.csv and return a DataFrame with the columns we need.

    Row order is preserved so that strided sharding (iloc[shard_id::num_shards])
    is deterministic across processes.
    """
    df = pd.read_csv(input_csv, low_memory=False)
    missing = [c for c in INPUT_KEY_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"input {input_csv} missing columns: {missing}")
    return df


def make_variant(row):
    """pooled_labeled.csv row -> genome.Variant (1-based, forward strand)."""
    chrom = str(row["chrom"])
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return genome.Variant(
        chromosome=chrom,
        position=int(row["pos"]),
        reference_bases=str(row["ref"]),
        alternate_bases=str(row["alt"]),
        name=str(row["accession"]),
    )


def make_interval(variant, seq_length_key=DEFAULT_SEQ_LENGTH):
    """Variant-centred interval resized to a supported model length."""
    return variant.reference_interval.resize(SEQUENCE_LENGTHS[seq_length_key])


# Scoring (one variant), with a small retry for transient gRPC errors

def score_one(client, variant, scorers, seq_length_key=DEFAULT_SEQ_LENGTH,
              retries=3, backoff=2.0):
    """Score a single variant and return AlphaGenome's tidy long DataFrame.

    Raises the last exception if all retries fail (caller decides what to do).
    """
    interval = make_interval(variant, seq_length_key)
    last_err = None
    for attempt in range(retries):
        try:
            scores = client.score_variant(interval, variant, variant_scorers=scorers)
            return variant_scorers.tidy_scores(scores)
        except Exception as e:  # noqa: BLE001 - includes grpc.RpcError
            last_err = e
            if attempt < retries - 1:
                time.sleep(backoff * (attempt + 1))
    raise last_err


# Reduction: tidy long -> one wide row per variant

def _sanitize(name):
    out = []
    for ch in str(name):
        out.append(ch if (ch.isalnum()) else "_")
    s = "".join(out)
    while "__" in s:
        s = s.replace("__", "_")
    return s.strip("_")


def _summ(values, prefix, out):
    """Signed + magnitude summaries of a score vector into `out` dict."""
    v = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy()
    out[f"{prefix}_n"] = int(v.size)
    if v.size == 0:
        out[f"{prefix}_mean"] = np.nan
        out[f"{prefix}_maxabs"] = np.nan
        out[f"{prefix}_signed_at_maxabs"] = np.nan
        return
    out[f"{prefix}_mean"] = float(np.mean(v))
    i = int(np.argmax(np.abs(v)))
    out[f"{prefix}_maxabs"] = float(np.abs(v[i]))
    out[f"{prefix}_signed_at_maxabs"] = float(v[i])


def _gene_subset(df, gene):
    """Prefer rows for the target gene; fall back to all rows if that scorer
    produced no gene-labelled rows (e.g. non-gene-centric scorers)."""
    if "gene_name" in df.columns:
        sub = df[df["gene_name"].astype(str).str.upper() == str(gene).upper()]
        if len(sub):
            return sub
    return df


def reduce_to_wide(tidy, gene, accession):
    """Collapse one variant's tidy long table into a single wide-row dict.

    `tidy` is the output of tidy_scores for ONE variant (all six scorers).
    `gene` filters to the assayed gene (a 1 MB window overlaps neighbours).
    Version note: v1 reduction, validated against pilot output before the full run.
    """
    out = {"accession": accession, "ag_gene": gene}
    if tidy is None or len(tidy) == 0:
        out["ag_status"] = "empty"
        return out
    out["ag_status"] = "ok"
    out["ag_n_rows_total"] = int(len(tidy))

    vs_col = tidy["variant_scorer"].astype(str) if "variant_scorer" in tidy else pd.Series([""] * len(tidy))
    ot_col = tidy["output_type"].astype(str) if "output_type" in tidy else pd.Series([""] * len(tidy))
    raw = "raw_score"
    quant = "quantile_score" if "quantile_score" in tidy.columns else None

    # ---- Expression: GeneMaskLFCScorer, RNA_SEQ (signed log2FC) ----
    lfc = tidy[vs_col.str.contains("LFC", na=False) & (ot_col == "RNA_SEQ")]
    lfc = _gene_subset(lfc, gene)
    _summ(lfc[raw] if len(lfc) else [], "expr_lfc", out)
    if quant and len(lfc):
        out["expr_lfc_quantile_maxabs"] = float(pd.to_numeric(lfc[quant], errors="coerce").abs().max())

    # tissue-resolved expression: one column per GTEx tissue (mean over its tracks).
    # Only RNA-seq tracks that are actually GTEx-derived carry a gtex_tissue; the
    # rest (e.g. ENCODE cell lines) leave it blank/NaN and must be excluded, else
    # they collapse into one spurious empty-named tissue and inflate n_gtex_tissues.
    n_gtex = 0
    if len(lfc) and "gtex_tissue" in lfc.columns:
        tis = lfc["gtex_tissue"].astype(str).str.strip()
        g = lfc[tis.ne("") & tis.str.lower().ne("nan")]
        if len(g):
            per_tissue = (
                g.assign(_r=pd.to_numeric(g[raw], errors="coerce"))
                 .groupby("gtex_tissue")["_r"].mean()
            )
            n_gtex = int(per_tissue.notna().sum())
            for tissue, val in per_tissue.items():
                out[f"expr_lfc_gtex__{_sanitize(tissue)}"] = float(val)
    out["n_gtex_tissues"] = n_gtex
    out["has_gtex_track"] = bool(n_gtex > 0)

    # ---- Expression "active" magnitude (GeneMaskActiveScorer, RNA_SEQ) ----
    act = tidy[vs_col.str.contains("Active", na=False) & (ot_col == "RNA_SEQ")]
    act = _gene_subset(act, gene)
    _summ(act[raw] if len(act) else [], "expr_active", out)

    # ---- Splicing ----
    ss = _gene_subset(tidy[ot_col == "SPLICE_SITES"], gene)
    _summ(ss[raw] if len(ss) else [], "splice_sites", out)
    su = _gene_subset(tidy[ot_col == "SPLICE_SITE_USAGE"], gene)
    _summ(su[raw] if len(su) else [], "splice_usage", out)
    sj = _gene_subset(tidy[ot_col == "SPLICE_JUNCTIONS"], gene)
    _summ(sj[raw] if len(sj) else [], "splice_junctions", out)

    # ---- Polyadenylation (UTR/abundance branch) ----
    pa = tidy[vs_col.str.contains("Polyadenylation", na=False)]
    pa = _gene_subset(pa, gene)
    _summ(pa[raw] if len(pa) else [], "polya", out)

    return out
