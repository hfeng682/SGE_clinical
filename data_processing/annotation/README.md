# annotation — raw SGE score files → annotated variants

Reproduces the transform from **raw MaveDB two-score SGE score files** to the
**annotated files** used downstream (genomic coordinates, harmonized consequence
labels, protein-visible split, distance-to-junction bins).

Pure standard library + pandas. Genomic annotation is fetched live from the
Ensembl VEP REST API, so the machine must have network access to
`rest.ensembl.org`.

## Layout

| module | responsibility |
|---|---|
| `config.py` | per-gene settings (transcript, strand, chrom, RNA column) + bin edges |
| `vep.py` | batch Ensembl VEP: transcript-HGVS → genomic pos, consequence, amino acids |
| `normalize.py` | forward-genomic ref/alt (strand flip) + VCF left-alignment of indels |
| `consequence.py` | coarse consequence vocabulary + route_class + protein_visible |
| `distance.py` | HGVS-offset parsing, exon-boundary model, distance to nearest junction |
| `pipeline.py` | orchestrates the five stages, writes the annotated CSV |

## Usage

```bash
# one gene
python -m annotation.pipeline BARD1 raw/BARD1.csv out/BARD1_annotated.csv

# all five configured genes (expects raw/<GENE>.csv)
python -m annotation.pipeline --all raw/ out/
```

```python
from annotation.pipeline import transform
df = transform("VHL", "raw/VHL.csv", "out/VHL_annotated.csv")
```

## Pipeline stages

1. **VEP annotation** — each raw `hgvs_nt` is rewritten onto the configured
   transcript and POSTed to VEP in batches of 200. Genomic chrom/pos,
   consequence terms, and amino-acid change are taken from the **target
   transcript only**. (BRCA1 note: the raw file uses `NM_007294.3`, which is
   retired at Ensembl; `config.py` bumps it to `NM_007294.4`, identical coding
   coordinates.)
2. **Forward-genomic alleles + indel left-alignment** — VEP returns alleles on
   the *transcript* strand, so on minus-strand genes (BARD1/PALB2/RAD51D/BRCA1)
   `ref`/`alt` are reverse-complemented to forward-genomic. Indels are then
   left-aligned to VCF convention (matching ClinVar/gnomAD). VHL is plus-strand
   and substitution-only.
3. **Consequence harmonization** — `coarse_consequence` collapses VEP's fine
   terms onto a 10-term vocabulary shared by all genes; `route_class` and
   `protein_visible` implement the protein-visible / protein-blind split.
4. **Distance to junction** — intronic distance is the HGVS offset itself;
   exonic distance is measured to the nearest exon edge, where edges are
   reconstructed from the intron offsets present in the (saturation) data.
   `pipeline.py` prints a donor→acceptor(+1) pairing count as a self-check;
   unpaired boundaries indicate screen coverage gaps (expected for
   non-contiguous exon designs, e.g. BRCA1).
5. **Column ordering + write.**

## Output schema

Constant annotation columns:

```
accession, hgvs_nt, chrom, pos, ref, alt,
consequence, coarse_consequence, route_class, protein_visible,
amino_acids, protein_pos, most_severe,
<original per-screen score / RNA columns>,
dist_to_junction, junction_side, dist_bin
```

The per-screen score columns are passed through unchanged (e.g. BARD1 keeps
`score`, `standard_error`, GMM densities, `rna_score`; BRCA1 keeps its replicate
columns; VHL keeps its two RNA timepoints). The mRNA-abundance score differs by
file — see `rna_score_col` in `config.py`.

## Validation (this release)

The package output was checked against the in-session annotated files:

* row-local transforms (pos, forward ref/alt incl. left-aligned indels, coarse
  consequence, route_class, protein_visible) — **100%** on stratified samples of
  all five genes;
* distance-to-junction on a **full** RAD51D rerun — **5410/5410** exact for both
  `dist_to_junction` and `dist_bin`.

Distances are validated *within* this session; the coordinate annotation was
independently cross-checked against the older coordinate-bearing MaveDB releases
(BARD1/PALB2/RAD51D) during development (100% on genomic position and
amino-acid change for shared variants).

## Adding a new screen

Append one entry to `GENES` in `config.py` (transcript, refseq flag, chrom,
strand, RNA column). Nothing else needs to change — reference sequence for
indel alignment is fetched on demand.
