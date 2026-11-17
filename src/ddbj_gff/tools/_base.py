"""Shared CLI skeleton for the per-tool converters (python -m ddbj_gff.tools.<tool>).

A tool module is a ToolSpec: which reader, whether to run the ID-dedupe prepass,
which NormalizeConfig flags to set and which normalize passes to run, in order.
Output is a plain GFF3 (only ##gff-version is guaranteed); run normalize afterwards.
"""
from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass, field

from .. import parse
from ..gtf import parse_gtf
from ..io import open_text
from ..model import Directive
from ..normalize import passes as _passes
from ..normalize.config import NormalizeConfig
from ..normalize.normalize import _APPLIED
from ..normalize.passes import NormalizeContext
from ..normalize.report import NormalizationReport
from ..validate.vocab import load_vocab
from ..writer import write
from .prepass import _is_gtf, dedupe_ids_per_seqid
from .snap_reader import parse_snap

_READERS = {"gff3": parse, "gtf": parse_gtf, "snap": parse_snap}


@dataclass
class ToolSpec:
    name: str
    reader: str = "auto"                 # "gff3" | "gtf" | "snap" | "auto" (.gtf suffix -> gtf)
    dedupe_ids: bool = False
    config: dict = field(default_factory=dict)
    passes: tuple[str, ...] = ()


def _reader_kind(spec: ToolSpec, text: str) -> str:
    if spec.reader != "auto":
        return spec.reader
    return "gtf" if _is_gtf(text.split("\n")) else "gff3"      # content, not suffix (AUGUSTUS writes GTF as .gff)


def convert(spec: ToolSpec, text: str) -> tuple:
    """Apply the spec to GFF/GTF/SNAP text. Returns (doc, report, n_renamed)."""
    n_renamed = 0
    if spec.dedupe_ids:
        text, n_renamed = dedupe_ids_per_seqid(text)
    doc = _READERS[_reader_kind(spec, text)](text)
    if doc.gff_version is None:
        doc.directives.insert(0, Directive("##gff-version 3", "gff-version", "3"))
    cfg = NormalizeConfig(**spec.config)
    ctx = NormalizeContext(vocab=load_vocab(), seq_lengths=None, config=cfg)
    report = NormalizationReport()
    for name in spec.passes:
        for change in getattr(_passes, name)(doc, ctx):
            (report.applied if change.action in _APPLIED else report.unresolved).append(change)
    return doc, report, n_renamed


def _reader_summary(doc) -> str:
    """One line per diagnostic code the reader raised (rows it skipped or synthesized)."""
    by_code: dict[str, list] = {}
    for d in doc.diagnostics:
        by_code.setdefault(d.code, []).append(d)
    return "".join(f"reader: {code} {len(ds)} ({ds[0].severity.value}): {ds[0].message}\n"
                   for code, ds in by_code.items())


def run_tool(spec: ToolSpec, argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog=f"gff_from_tool {spec.name}",
                                 description=f"Convert {spec.name} output to a standard gene/mRNA/CDS GFF3")
    ap.add_argument("input")
    ap.add_argument("--out")
    ap.add_argument("--report")
    args = ap.parse_args(argv)
    if not os.path.exists(args.input):
        sys.stderr.write(f"error: input file not found: {args.input}\n")
        return 2
    with open_text(args.input, errors="replace") as fh:
        text = fh.read()
    doc, report, n_renamed = convert(spec, text)

    out_text = write(doc)
    if args.out:
        with open(args.out, "w", encoding="ascii", errors="replace") as fh:
            fh.write(out_text)
    else:
        sys.stdout.write(out_text)
    report_text = report.render() + _reader_summary(doc)
    if spec.dedupe_ids:
        report_text += f"prepass: renamed {n_renamed} id(s) that collided across seqids\n"
    if args.report:
        with open(args.report, "w", encoding="utf-8") as fh:
            fh.write(report_text)
    else:
        sys.stderr.write(report_text)
    return 0
