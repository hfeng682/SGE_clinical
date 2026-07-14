"""
Consequence harmonization.

Two derived columns from VEP's fine-grained consequence terms:

  coarse_consequence : collapses VEP terms onto the 10-term vocabulary used by
                       the MaveDB '-a-1' annotated releases, so all five genes
                       share one label set. Validated in-session to agree 100%
                       with the old files on shared variants.

  route_class + protein_visible : the project's protein-visible vs protein-blind
                       split. protein_visible variants are those whose effect can
                       act through the encoded protein (missense/nonsense/
                       synonymous/indel/start-stop-loss); protein-blind variants
                       (intronic, splice_region, UTR) can only act through RNA.
"""
from __future__ import annotations

PROTEIN_VISIBLE = {"missense", "nonsense", "synonymous", "inframe/indel", "start/stop_lost"}


def coarse_consequence(cons: str):
    """VEP fine consequence (';'-joined) -> single coarse label (old-file vocabulary)."""
    if not isinstance(cons, str):
        return None
    s = set(cons.split(";"))
    has = lambda *t: any(x in s for x in t)
    if has("splice_acceptor_variant", "splice_donor_variant"): return "splice_site_variant"
    if has("stop_gained"):                                     return "stop_gained"
    if has("start_lost"):                                      return "start_lost"
    if has("stop_lost"):                                       return "stop_lost"
    if has("missense_variant"):                                return "missense_variant"
    if has("frameshift_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant"):
        return "inframe_indel"
    if any("splice" in x for x in s):                          return "splicing_variant"
    if has("synonymous_variant", "stop_retained_variant", "start_retained_variant"):
        return "synonymous_variant"
    if has("3_prime_UTR_variant", "5_prime_UTR_variant"):      return "UTR_variant"
    if has("intron_variant"):                                  return "intron_variant"
    return sorted(s)[0] if s else None


def route_class(cons: str):
    """Finer project stratum used for the protein/RNA route logic."""
    if not isinstance(cons, str):
        return "unknown"
    s = set(cons.split(";"))
    if "missense_variant" in s:                                return "missense"
    if "stop_gained" in s:                                     return "nonsense"
    if s & {"synonymous_variant", "stop_retained_variant", "start_retained_variant"}:
        return "synonymous"
    if any("splice" in x for x in s):                          return "splice_region"
    if s & {"frameshift_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant"}:
        return "inframe/indel"
    if "intron_variant" in s:                                  return "intronic"
    if "5_prime_UTR_variant" in s:                             return "5UTR"
    if "3_prime_UTR_variant" in s:                             return "3UTR"
    if s & {"start_lost", "stop_lost"}:                        return "start/stop_lost"
    return sorted(s)[0] if s else "unknown"


def is_protein_visible(rclass: str) -> bool:
    return rclass in PROTEIN_VISIBLE
