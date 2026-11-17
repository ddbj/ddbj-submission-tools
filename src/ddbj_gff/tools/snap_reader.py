"""Reader for SNAP output as written by zff2gff3.pl: 8-column rows
(seqid source CDS start end score strand Name=<gene>) — one row per exon, grouped
into gene -> mRNA -> CDS by Name. Phases are left unset (pass_compute_phase)."""
from __future__ import annotations

from ..attributes import parse_attributes
from ..errors import Diagnostic, Severity
from ..model import Feature, GffDocument, Span
from ..parser import _resolve_graph


def parse_snap(text: str) -> GffDocument:
    doc = GffDocument()
    groups: dict[str, list[tuple[str, Span]]] = {}
    for line_no, line in enumerate(text.splitlines(), start=1):
        if not line or line.startswith("#"):
            continue
        cols = line.rstrip("\r\n").split("\t")
        if len(cols) == 8:
            attr_s = cols[7]
        elif len(cols) == 9:
            attr_s = cols[8]
        else:
            doc.diagnostics.append(Diagnostic(Severity.ERROR, line_no, "col-count",
                                              f"expected 8 or 9 columns, got {len(cols)}"))
            continue
        seqid, source, ftype, start_s, end_s, score_s, strand = cols[:7]
        name = parse_attributes(attr_s).get("Name", [None])[0]
        if ftype != "CDS" or name is None:
            doc.diagnostics.append(Diagnostic(Severity.WARNING, line_no, "snap-row-skipped",
                                              f"{ftype} row without Name; skipped"))
            continue
        score = None if score_s == "." else float(score_s)
        groups.setdefault(name, []).append((source, Span(seqid, int(start_s), int(end_s), strand, None, score)))

    for name, rows in groups.items():
        source = rows[0][0]
        spans = [s for _, s in rows]
        lo, hi = min(s.start for s in spans), max(s.end for s in spans)
        hull = (spans[0].seqid, lo, hi, spans[0].strand)
        tid, cid = f"{name}.t1", f"{name}.t1.cds"
        gene = Feature(name, source, "gene", [Span(*hull)], {"ID": [name]}, [])
        mrna = Feature(tid, source, "mRNA", [Span(*hull)], {"ID": [tid], "Parent": [name]}, [name])
        cds = Feature(cid, source, "CDS", spans, {"ID": [cid], "Parent": [tid]}, [tid])
        for f in (gene, mrna, cds):
            doc.features.append(f)
            doc.feature_index[f.id] = f
    _resolve_graph(doc)
    return doc
