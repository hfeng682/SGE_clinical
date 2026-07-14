"""
Coordinate / allele normalization.

Two jobs:
  1. forward-genomic ref/alt  -- VEP reports alleles on the transcript strand;
     on minus-strand genes we reverse-complement to forward-genomic so the
     variants join against ClinVar/gnomAD and against each other on chrom:pos:ref:alt.
  2. left-normalized indels   -- standard VCF left-alignment (the convention
     ClinVar/gnomAD use). SNV-only screens (BRCA1, VHL) never hit this path.

Reference sequence for indel left-alignment is fetched once per gene region from
the Ensembl REST /sequence endpoint and cached.
"""
from __future__ import annotations
import json, urllib.request, urllib.parse
from .config import ENSEMBL_REST

_COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "-": "-"}


def revcomp(s: str) -> str:
    return "".join(_COMP.get(b, b) for b in str(s)[::-1])


def forward_alleles(allele_string: str, strand: int):
    """VEP allele_string ('REF/ALT', transcript strand) -> (ref, alt) forward-genomic."""
    if not isinstance(allele_string, str) or "/" not in allele_string:
        return (None, None)
    parts = allele_string.split("/")
    ref_t, alt_t = parts[0], parts[-1]
    if strand == -1:
        return (revcomp(ref_t), revcomp(alt_t))
    return (ref_t, alt_t)


class RefSeqCache:
    """Fetch and cache a padded reference window covering a gene's variants."""
    def __init__(self):
        self._regions = {}   # gene -> (lo, hi, seq)

    def load(self, gene, chrom, pos_min, pos_max, pad=50):
        lo, hi = int(pos_min) - pad, int(pos_max) + pad
        url = f"{ENSEMBL_REST}/sequence/region/human/{chrom}:{lo}..{hi}:1?content-type=application/json"
        with urllib.request.urlopen(url, timeout=60) as r:
            seq = json.load(r)["seq"].upper()
        self._regions[gene] = (lo, hi, seq)

    def base(self, gene, pos):
        lo, hi, seq = self._regions[gene]
        return seq[int(pos) - lo]

    def has(self, gene):
        return gene in self._regions


def left_align(cache: RefSeqCache, gene, pos, ref, alt):
    """
    Canonical bcftools-style VCF left-alignment.
    ref/alt may be '' or '-' for the empty side of an indel. Returns (pos, ref, alt).
    SNVs are returned unchanged.
    """
    ref = "" if ref in ("-", None) else str(ref)
    alt = "" if alt in ("-", None) else str(alt)
    pos = int(pos)
    if len(ref) == 1 and len(alt) == 1:
        return pos, ref, alt
    while True:
        changed = False
        if ref and alt and ref[-1] == alt[-1]:          # trim shared suffix
            ref, alt = ref[:-1], alt[:-1]; changed = True
        if len(ref) == 0 or len(alt) == 0:              # controls + roll left
            b = cache.base(gene, pos - 1)
            ref, alt, pos = b + ref, b + alt, pos - 1; changed = True
        if not changed:
            break
    return pos, ref, alt
