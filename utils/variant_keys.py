"""
variant_keys.py -- the single variant identity used to join every table.

A variant is keyed on (gene, chrom, pos, ref, alt) on GRCh38, VCF-left-aligned.
All merges in the repo go through `key()` so the join is defined in exactly one
place.
"""
from __future__ import annotations
import pandas as pd

KEY = ["gene", "chrom", "pos", "ref", "alt"]

def key(df: pd.DataFrame) -> pd.Series:
    """Row-wise '|'-joined variant key over the frozen KEY columns."""
    return df[KEY].astype(str).agg("|".join, axis=1)
