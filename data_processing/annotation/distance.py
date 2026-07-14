"""
Distance to nearest splice junction, derived from the HGVS c. notation.

Design:
  * Intronic variants carry the offset explicitly (c.1903+62 = 62 bp from that
    donor; c.302-15 = 15 bp from that acceptor). The offset IS the distance.
  * Exonic variants have no offset, so distance is measured to the nearest exon
    edge. Because these are saturation screens with dense intronic coverage,
    every intron offset pins an exon boundary in a linear cDNA coordinate
    (5'UTR -> CDS -> 3'UTR); we reconstruct those boundaries from the data
    itself and measure each exonic variant to the closest one.

The exon model is self-checking: each donor boundary d should pair with an
acceptor boundary at d+1 (see pipeline.py, which reports the pairing count).
"""
from __future__ import annotations
import re
import numpy as np


def parse_c(hgvs_nt: str):
    """
    Parse the c. portion of an hgvs_nt string.
    Returns dict(seg, base, offset): seg in {'5utr','cds','3utr'}, base is the
    coding coordinate (negative for 5'UTR), offset is the intron offset (0 exonic).
    """
    if ":c." not in str(hgvs_nt):
        return dict(seg=None, base=None, offset=None)
    v = hgvs_nt.split(":c.")[1]
    m = re.match(r"([*-]?)(\d+)([+-]\d+)?", v)
    if not m:
        return dict(seg=None, base=None, offset=None)
    pre, num, off = m.group(1), int(m.group(2)), m.group(3)
    offset = int(off) if off else 0
    if pre == "*":
        seg = "3utr"
    elif pre == "-":
        seg, num = "5utr", -num
    else:
        seg = "cds"
    return dict(seg=seg, base=num, offset=offset)


def linear_coord(seg, base, cds_len):
    """Place a variant on a single linear cDNA axis: 5'UTR(<0) -> CDS -> 3'UTR(>cds_len)."""
    if seg == "3utr":
        return cds_len + int(base)
    return int(base)             # cds positive, 5utr already negative


def build_exon_model(parsed_rows):
    """
    parsed_rows: list of dicts from parse_c (one per variant).
    Returns dict(cds_len, donors, acceptors, boundaries) in linear cDNA coordinates.
      donors    = exon-end edges  (positions with a +offset intronic neighbour)
      acceptors = exon-start edges (positions with a -offset intronic neighbour)
    """
    cds_bases = [r["base"] for r in parsed_rows if r["seg"] == "cds" and r["offset"] == 0]
    cds_len = max(cds_bases) if cds_bases else 0
    donors, acceptors = set(), set()
    for r in parsed_rows:
        if r["seg"] is None or r["offset"] == 0:
            continue
        lin = linear_coord(r["seg"], r["base"], cds_len)
        if r["offset"] > 0:
            donors.add(lin)
        else:
            acceptors.add(lin)
    boundaries = np.array(sorted(donors | acceptors), dtype=float)
    return dict(cds_len=cds_len, donors=donors, acceptors=acceptors, boundaries=boundaries)


def distance_and_side(parsed_row, model):
    """Return (distance_bp, side) for one variant. side in {'donor','acceptor', None}."""
    seg, base, off = parsed_row["seg"], parsed_row["base"], parsed_row["offset"]
    if seg is None:
        return (np.nan, None)
    if off != 0:                                        # intronic: offset is the distance
        return (abs(off), "donor" if off > 0 else "acceptor")
    # exonic: nearest boundary in linear coordinate
    b = model["boundaries"]
    if len(b) == 0:
        return (np.nan, None)
    lin = linear_coord(seg, base, model["cds_len"])
    j = int(np.argmin(np.abs(b - lin)))
    edge = b[j]
    side = "donor" if edge in model["donors"] else "acceptor"
    return (abs(edge - lin), side)


def pairing_report(model):
    """Self-consistency: how many donor edges d pair with an acceptor at d+1."""
    acc = model["acceptors"]
    paired = sum(1 for d in model["donors"] if (d + 1) in acc)
    return paired, len(model["donors"])
