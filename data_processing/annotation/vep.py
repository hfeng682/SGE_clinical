"""
Ensembl VEP REST client: transcript-HGVS -> genomic position + consequence.

The raw MaveDB files give each variant only as a transcript-level HGVS string
(e.g. 'NM_000465.4:c.1670G>A'). This module posts those to the VEP 'hgvs'
endpoint in batches and extracts, for the target transcript only, the fields the
pipeline needs: genomic chrom/pos, VEP allele_string, consequence terms, and
amino-acid change.

Only the Python standard library is used (urllib), so the module has no
third-party dependencies.
"""
from __future__ import annotations
import json, time, urllib.request, urllib.parse, urllib.error
from .config import ENSEMBL_REST

BATCH = 200          # VEP hgvs POST accepts up to 200 notations
SLEEP = 0.34         # ~3 requests/s, within Ensembl's rate limit


def _post(hgvs_list, refseq, timeout=120, retries=5):
    """POST one batch of HGVS notations; return list of VEP records."""
    url = f"{ENSEMBL_REST}/vep/human/hgvs?content-type=application/json" + ("&refseq=1" if refseq else "")
    data = json.dumps({"hgvs_notations": list(hgvs_list)}).encode()
    req = urllib.request.Request(
        url, data=data,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as e:
            if e.code == 429 or e.code >= 500:       # rate-limited / transient
                time.sleep(2 + 2 * attempt); continue
            raise                                    # 400 etc. -> surface
        except Exception:
            time.sleep(2 + attempt)
    return []


def _parse(record, transcript):
    """Pull target-transcript fields out of one VEP record."""
    txbase = transcript.split(".")[0]
    out = dict(
        query=record.get("input"),
        chrom=record.get("seq_region_name"),
        pos=record.get("start"),
        allele_string=record.get("allele_string"),
        most_severe=record.get("most_severe_consequence"),
        consequence=None, amino_acids=None, protein_pos=None,
    )
    for c in record.get("transcript_consequences", []):
        if (c.get("transcript_id") or "").startswith(txbase):
            out.update(
                consequence=";".join(c.get("consequence_terms", [])),
                amino_acids=c.get("amino_acids"),
                protein_pos=c.get("protein_start"),
            )
            break
    return out


def build_query(transcript, hgvs_nt):
    """Rewrite a raw hgvs_nt onto the configured transcript (handles BRCA1 .3->.4)."""
    return transcript + ":" + hgvs_nt.split(":", 1)[1]


def annotate(hgvs_nt_series, transcript, refseq, progress=False):
    """
    Annotate a series/list of raw hgvs_nt strings.
    Returns a list of parsed dicts aligned 1:1 with the input order.
    """
    queries = [build_query(transcript, h) for h in hgvs_nt_series]
    parsed = {}
    for i in range(0, len(queries), BATCH):
        chunk = queries[i:i + BATCH]
        for rec in _post(chunk, refseq):
            parsed[rec.get("input")] = _parse(rec, transcript)
        if progress:
            print(f"  VEP {min(i + BATCH, len(queries))}/{len(queries)}")
        time.sleep(SLEEP)
    # align to input order; queries are unique per row within a gene
    return [parsed.get(q, {"query": q}) for q in queries]
