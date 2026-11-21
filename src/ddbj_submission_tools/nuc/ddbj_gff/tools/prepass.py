"""Text-level prepass run before parsing: rename IDs that collide across seqids.

Per-chromosome runs concatenated into one file (AUGUSTUS `g1` on every chromosome)
would be merged by the parser into one bogus discontinuous feature. Rewriting such
IDs to `<seqid>.<ID>` (and the Parent= / gene_id / transcript_id that reference them
on the same seqid) keeps every model distinct. IDs repeated only within one seqid
(discontinuous CDS rows) are untouched; so are comment lines and the ##FASTA section.
"""
from __future__ import annotations

import re

_GFF_ID = re.compile(r"(?<![\w.])(ID|Parent)=([^;\n]+)")
_GTF_ID = re.compile(r'(gene_id|transcript_id) "([^"]+)"')


def _is_gtf(lines: list[str]) -> bool:
    for line in lines:
        if line.startswith("##FASTA"):
            break
        if not line or line.startswith("#"):
            continue
        if 'gene_id "' in line or 'transcript_id "' in line:
            return True
        if "ID=" in line or "Parent=" in line:
            return False
    return False


def _ids_in(line: str, gtf: bool) -> list[str]:
    cols = line.split("\t")
    if len(cols) < 9:
        return []
    col9 = cols[8]
    if gtf:
        if '"' not in col9 and cols[2] in ("gene", "transcript"):
            return [col9.strip()]
        return [v for _, v in _GTF_ID.findall(col9)]
    return [v for k, vs in _GFF_ID.findall(col9) for v in vs.split(",") if k == "ID"] + \
           [v for k, vs in _GFF_ID.findall(col9) for v in vs.split(",") if k == "Parent"]


def dedupe_ids_per_seqid(text: str) -> tuple[str, int]:
    lines = text.split("\n")
    gtf = _is_gtf(lines)
    seqids_of: dict[str, set[str]] = {}
    body_end = len(lines)
    for i, line in enumerate(lines):
        if line.startswith("##FASTA"):
            body_end = i
            break
        if not line or line.startswith("#"):
            continue
        seqid = line.split("\t", 1)[0]
        for v in _ids_in(line, gtf):
            seqids_of.setdefault(v, set()).add(seqid)
    dup = {v for v, s in seqids_of.items() if len(s) > 1}
    if not dup:
        return text, 0

    def gff_sub(m, seqid):
        key, vals = m.group(1), m.group(2)
        new = ",".join(f"{seqid}.{v}" if v in dup else v for v in vals.split(","))
        return f"{key}={new}"

    def gtf_sub(m, seqid):
        key, v = m.group(1), m.group(2)
        return f'{key} "{seqid}.{v}"' if v in dup else m.group(0)

    out = list(lines)
    for i in range(body_end):
        line = lines[i]
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        seqid = cols[0]
        if gtf and '"' not in cols[8] and cols[2] in ("gene", "transcript"):
            v = cols[8].strip()
            if v in dup:
                cols[8] = f"{seqid}.{v}"
        elif gtf:
            cols[8] = _GTF_ID.sub(lambda m: gtf_sub(m, seqid), cols[8])
        else:
            cols[8] = _GFF_ID.sub(lambda m: gff_sub(m, seqid), cols[8])
        out[i] = "\t".join(cols)
    return "\n".join(out), len(dup)
