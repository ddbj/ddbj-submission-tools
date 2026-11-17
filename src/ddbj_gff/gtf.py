"""GTF (2.2) reader: builds the same GffDocument as parse() would for an
equivalent GFF3 (gene -> transcript -> exon/CDS/...), from gene_id/transcript_id."""
from __future__ import annotations

from .errors import Diagnostic, Severity
from .model import Feature, GffDocument, Span
from .parser import _resolve_graph

_PARENT_ROWS = {"gene", "transcript", "mRNA"}


def parse_gtf_attributes(col9: str) -> dict[str, list[str]]:
    attrs: dict[str, list[str]] = {}
    for item in col9.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        key, _, val = item.partition(" ")
        val = val.strip().strip('"')
        attrs.setdefault(key, []).append(val)
    return attrs


def parse_gtf(text: str) -> GffDocument:
    doc = GffDocument()
    # per-transcript bookkeeping (insertion-ordered)
    genes: dict[str, Feature] = {}
    transcripts: dict[str, Feature] = {}
    tx_gene: dict[str, str] = {}
    cds: dict[str, Feature] = {}
    counters: dict[tuple[str, str], int] = {}
    children: list[Feature] = []
    n_syn_gene = n_syn_tx = 0

    for line_no, line in enumerate(text.splitlines(), start=1):
        if not line or line.startswith("#"):
            continue
        cols = line.rstrip("\r\n").split("\t")
        if len(cols) != 9:
            doc.diagnostics.append(Diagnostic(Severity.ERROR, line_no, "col-count",
                                              f"expected 9 columns, got {len(cols)}"))
            continue
        seqid, source, ftype, start_s, end_s, score_s, strand, phase_s, col9 = cols
        span = Span(seqid, int(start_s), int(end_s), strand,
                    None if phase_s == "." else int(phase_s),
                    None if score_s == "." else float(score_s))
        if '"' not in col9 and ftype in ("gene", "transcript"):
            # AUGUSTUS/BRAKER GTF: bare token id
            attrs = {"gene_id": [col9.strip()]} if ftype == "gene" else {"transcript_id": [col9.strip()]}
            if ftype == "transcript":
                attrs["gene_id"] = [col9.strip().rsplit(".", 1)[0]]
        else:
            attrs = parse_gtf_attributes(col9)
        gid = attrs.pop("gene_id", [None])[0]
        tid = attrs.pop("transcript_id", [None])[0]

        if ftype == "gene":
            if gid is not None and gid not in genes:
                genes[gid] = Feature(gid, source, "gene", [span], {"ID": [gid], **attrs}, [])
            continue
        if ftype in _PARENT_ROWS:
            if tid is not None and tid not in transcripts:
                transcripts[tid] = Feature(tid, source, ftype, [span],
                                           {"ID": [tid], "Parent": [gid], **attrs}, [gid])
                tx_gene[tid] = gid
            continue
        if tid is None:
            doc.diagnostics.append(Diagnostic(Severity.WARNING, line_no, "gtf-no-transcript-id",
                                              f"{ftype} row without transcript_id; skipped"))
            continue
        tx_gene.setdefault(tid, gid)
        if ftype == "CDS":
            if tid not in cds:
                fid = f"{tid}.cds"
                cds[tid] = Feature(fid, source, "CDS", [], {"ID": [fid], "Parent": [tid], **attrs}, [tid])
                children.append(cds[tid])
            cds[tid].spans.append(span)
            continue
        n = counters[(tid, ftype)] = counters.get((tid, ftype), 0) + 1
        fid = f"{tid}.{ftype}{n}"
        children.append(Feature(fid, source, ftype, [span], {"ID": [fid], "Parent": [tid], **attrs}, [tid]))

    # synthesize missing transcript / gene rows from their children
    by_tx: dict[str, list[Feature]] = {}
    for c in children:
        by_tx.setdefault(c.parent_ids[0], []).append(c)
    for tid, kids in by_tx.items():
        if tid not in transcripts:
            gid = tx_gene.get(tid)
            transcripts[tid] = Feature(tid, kids[0].source, "transcript", [_hull(kids)],
                                       {"ID": [tid], "Parent": [gid]}, [gid])
            n_syn_tx += 1
    by_gene: dict[str, list[Feature]] = {}
    for tid, t in transcripts.items():
        by_gene.setdefault(t.parent_ids[0], []).append(t)
    for gid, txs in by_gene.items():
        if gid not in genes:
            genes[gid] = Feature(gid, txs[0].source, "gene", [_hull(txs)], {"ID": [gid]}, [])
            n_syn_gene += 1
    if n_syn_gene or n_syn_tx:
        doc.diagnostics.append(Diagnostic(Severity.INFO, None, "gtf-synthesized",
                                          f"{n_syn_gene} gene and {n_syn_tx} transcript rows synthesized"))

    for f in [*genes.values(), *transcripts.values(), *children]:
        doc.features.append(f)
        doc.feature_index[f.id] = f
    _resolve_graph(doc)
    return doc


def _hull(feats: list[Feature]) -> Span:
    spans = [s for f in feats for s in f.spans]
    return Span(spans[0].seqid, min(s.start for s in spans), max(s.end for s in spans), spans[0].strand)
