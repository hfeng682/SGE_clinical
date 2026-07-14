"""
Per-gene configuration for the two-score SGE annotation pipeline.

Each entry describes one MaveDB score set:
  transcript   : the reference transcript the pipeline annotates against.
                 NOTE the BRCA1 bump NM_007294.3 -> NM_007294.4: the .3 version
                 is retired at Ensembl and returns HTTP 400 from VEP; .4 maps
                 to identical coding coordinates.
  refseq       : True for RefSeq (NM_) transcripts  -> VEP called with refseq=1
                 False for Ensembl (ENST) transcripts.
  chrom        : GRCh38 chromosome (no 'chr' prefix; forward strand).
  strand       : gene strand (+1 / -1). VEP reports alleles on the TRANSCRIPT
                 strand, so on -1 genes the forward-genomic ref/alt are the
                 reverse complement of VEP's allele_string.
  rna_score_col: the column in the raw file that holds the mRNA-abundance score.
                 Canonical choice per screen; edit here if you prefer another.

These five genes are the datasets validated in-session. Adding a new screen is a
matter of appending one GENES entry (and, if it is on a new chromosome region,
nothing else — reference sequence is fetched on demand).
"""

GENES = {
    "BARD1":  dict(transcript="NM_000465.4",       refseq=True,  chrom="2",  strand=-1, rna_score_col="rna_score"),
    "PALB2":  dict(transcript="NM_024675.4",       refseq=True,  chrom="16", strand=-1, rna_score_col="rna_score"),
    "RAD51D": dict(transcript="NM_002878.4",       refseq=True,  chrom="17", strand=-1, rna_score_col="rna_score"),
    "BRCA1":  dict(transcript="NM_007294.4",       refseq=True,  chrom="17", strand=-1, rna_score_col="score_rna"),
    "VHL":    dict(transcript="ENST00000256474.3", refseq=False, chrom="3",  strand=1,  rna_score_col="rna_score_d20"),
}

ENSEMBL_REST = "https://rest.ensembl.org"
ASSEMBLY = "GRCh38"

# Distance-to-junction bin edges (bp). Exonic distance is measured to the nearest
# exon edge; intronic distance is the HGVS offset itself.
def distance_bin(offset_is_intronic: bool, dist: float) -> str:
    """Map (intronic?, distance_bp) -> categorical splice-proximity bin."""
    import math
    if dist is None or (isinstance(dist, float) and math.isnan(dist)):
        return "unknown"
    if not offset_is_intronic:                      # exonic
        return "exonic_splice_region(\u22643)" if dist <= 3 else "exonic_core(>3)"
    if dist <= 2:  return "splice_site(\u00b11-2)"  # essential dinucleotide
    if dist <= 8:  return "splice_region_intronic(3-8)"
    if dist <= 50: return "near_intronic(9-50)"
    return "deep_intronic(>50)"
