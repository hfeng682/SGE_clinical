"""
End-to-end transform: raw MaveDB score file -> annotated file.

Stages (per gene):
  1. VEP annotate raw hgvs_nt -> genomic pos, consequence, amino acids   (vep.py)
  2. forward-genomic ref/alt + VCF left-align indels                     (normalize.py)
  3. coarse_consequence + route_class + protein_visible                  (consequence.py)
  4. dist_to_junction + junction_side + dist_bin                         (distance.py)
  5. drop VEP columns that carry no data on this endpoint, order columns, write CSV

Output schema (score columns vary by screen; annotation columns are constant):
  accession, hgvs_nt, chrom, pos, ref, alt,
  consequence, coarse_consequence, route_class, protein_visible,
  amino_acids, protein_pos, most_severe,
  <original score/rna columns>,
  dist_to_junction, junction_side, dist_bin

Usage:
  python -m sge_pipeline.pipeline BARD1 raw/BARD1.csv out/BARD1_annotated.csv
  python -m sge_pipeline.pipeline --all raw/ out/
"""
from __future__ import annotations
import os, sys
import pandas as pd

from .config import GENES, distance_bin
from . import vep, normalize, consequence, distance

FRONT = ["accession", "hgvs_nt", "chrom", "pos", "ref", "alt",
         "consequence", "coarse_consequence", "route_class", "protein_visible",
         "amino_acids", "protein_pos", "most_severe"]
TAIL = ["dist_to_junction", "junction_side", "dist_bin"]


def transform(gene: str, raw_csv: str, out_csv: str, progress: bool = True) -> pd.DataFrame:
    if gene not in GENES:
        raise ValueError(f"unknown gene {gene!r}; known: {list(GENES)}")
    cfg = GENES[gene]
    df = pd.read_csv(raw_csv)
    if progress:
        print(f"[{gene}] {len(df)} variants from {raw_csv}")

    # ---- 1. VEP annotation -------------------------------------------------
    ann = vep.annotate(df["hgvs_nt"], cfg["transcript"], cfg["refseq"], progress=progress)
    A = pd.DataFrame(ann)
    df = df.reset_index(drop=True)
    df["chrom"] = A["chrom"]
    df["pos"] = pd.to_numeric(A["pos"], errors="coerce").astype("Int64")
    df["consequence"] = A["consequence"]
    df["amino_acids"] = A["amino_acids"]
    df["protein_pos"] = A["protein_pos"]
    df["most_severe"] = A["most_severe"]
    n_ann = df["consequence"].notna().sum()
    if progress:
        print(f"[{gene}] annotated on target transcript: {n_ann}/{len(df)}")

    # ---- 2. forward-genomic alleles + indel left-alignment -----------------
    fa = A["allele_string"].apply(lambda s: pd.Series(normalize.forward_alleles(s, cfg["strand"])))
    df["ref"], df["alt"] = fa[0], fa[1]
    is_indel = (df["ref"].astype(str).str.len() != 1) | (df["alt"].astype(str).str.len() != 1)
    if is_indel.any():
        cache = normalize.RefSeqCache()
        cache.load(gene, cfg["chrom"], df["pos"].min(), df["pos"].max())
        norm = [normalize.left_align(cache, gene, r.pos, r.ref, r.alt)
                for r in df[is_indel].itertuples()]
        nd = pd.DataFrame(norm, index=df[is_indel].index, columns=["pos", "ref", "alt"])
        df.loc[is_indel, ["pos", "ref", "alt"]] = nd
    df["pos"] = df["pos"].astype("Int64")

    # ---- 3. consequence harmonization --------------------------------------
    df["coarse_consequence"] = df["consequence"].map(consequence.coarse_consequence)
    df["route_class"] = df["consequence"].map(consequence.route_class)
    df["protein_visible"] = df["route_class"].map(consequence.is_protein_visible)

    # ---- 4. distance to junction -------------------------------------------
    parsed = [distance.parse_c(h) for h in df["hgvs_nt"]]
    model = distance.build_exon_model(parsed)
    paired, total = distance.pairing_report(model)
    if progress:
        print(f"[{gene}] exon model: {paired}/{total} donor->acceptor(+1) pairs "
              f"(unpaired = screen coverage gaps, expected for non-contiguous designs)")
    dd = [distance.distance_and_side(p, model) for p in parsed]
    df["dist_to_junction"] = [x[0] for x in dd]
    df["junction_side"] = [x[1] for x in dd]
    df["dist_bin"] = [distance_bin(p["offset"] != 0 if p["offset"] is not None else False, d)
                      for p, d in zip(parsed, df["dist_to_junction"])]

    # ---- 5. column ordering + write ----------------------------------------
    front = [c for c in FRONT if c in df.columns]
    tail = [c for c in TAIL if c in df.columns]
    middle = [c for c in df.columns if c not in front + tail]
    df = df[front + middle + tail]

    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    df.to_csv(out_csv, index=False)
    if progress:
        print(f"[{gene}] wrote {out_csv}  ({df.shape[0]} x {df.shape[1]})")
    return df


def _main(argv):
    if len(argv) >= 1 and argv[0] == "--all":
        raw_dir, out_dir = argv[1], argv[2]
        for gene in GENES:
            transform(gene, os.path.join(raw_dir, f"{gene}.csv"),
                      os.path.join(out_dir, f"{gene}_annotated.csv"))
    else:
        gene, raw_csv, out_csv = argv
        transform(gene, raw_csv, out_csv)


if __name__ == "__main__":
    _main(sys.argv[1:])
