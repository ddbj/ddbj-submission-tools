from __future__ import annotations

import argparse
import sys

from .. import parse
from ..writer import write
from .context import RepairContext
from .registry import list_operations, get_operation
from .report import candidates_to_json, render_candidates, render_changes
from .driver import run_detect, run_apply, load_sequences, apply_all_names


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="gff_repair",
                                 description="Apply modular GFF curation operations")
    ap.add_argument("--list", action="store_true", help="list available operations")
    ap.add_argument("--gff")
    ap.add_argument("--fasta")
    ap.add_argument("--transl-table", type=int, default=1, dest="transl_table")
    ap.add_argument("--inference-value", default="ab initio prediction",
                    dest="inference_value",
                    help="/inference value written by add-inference (verbatim, not validated)")
    ap.add_argument("--inference-types", default="mRNA", dest="inference_types",
                    help="comma-separated feature types for add-inference; '*' matches any")
    ap.add_argument("--inference-sources", default="", dest="inference_sources",
                    help="comma-separated GFF source values for add-inference; empty = any")
    ap.add_argument("--min-cds-length", type=int, default=90, dest="min_cds_length",
                    help="short-cds-drop-mrna threshold in bases (default 90, per ANN4300)")
    ap.add_argument("--detect", action="store_true", help="preview candidates; write nothing")
    ap.add_argument("--json", action="store_true", help="detect output as JSON")
    ap.add_argument("--only", help="comma-separated operation names (detect)")
    ap.add_argument("--apply", help="comma-separated operation names, or 'all'")
    ap.add_argument("--out")
    ap.add_argument("--report")
    args = ap.parse_args(argv)

    if args.list:
        for op in list_operations():
            seq = " (needs FASTA)" if op.requires_sequence else ""
            dest = " [destructive]" if op.destructive else ""
            sys.stdout.write(f"{op.name}{seq}{dest}: {op.summary}\n")
        return 0

    if not args.gff:
        ap.error("--gff is required unless --list")

    with open(args.gff, encoding="ascii", errors="replace") as fh:
        doc = parse(fh.read())
    sequences = load_sequences(args.fasta) if args.fasta else None
    ctx = RepairContext(sequences=sequences, transl_table=args.transl_table,
                        inference_value=args.inference_value,
                        inference_types=args.inference_types,
                        inference_sources=args.inference_sources,
                        min_cds_length=args.min_cds_length)

    if args.detect:
        names = args.only.split(",") if args.only else None
        try:
            ops = [get_operation(n) for n in names] if names else list_operations()
            missing_seq = sorted(op.name for op in ops
                                 if op.requires_sequence and ctx.sequences is None)
            if missing_seq:
                ap.error(f"operation(s) {', '.join(missing_seq)} require sequence data "
                        f"(pass --fasta)")
            cands = run_detect(doc, ctx, names)
        except KeyError as e:
            ap.error(str(e).strip('"'))
        sys.stdout.write(candidates_to_json(cands) + "\n" if args.json
                         else render_candidates(cands))
        return 0

    if args.apply:
        if args.apply == "all":
            names = apply_all_names()
        else:
            names = args.apply.split(",")
        try:
            ops = [get_operation(n) for n in names]
            missing_seq = sorted(op.name for op in ops
                                 if op.requires_sequence and ctx.sequences is None)
            if missing_seq:
                ap.error(f"operation(s) {', '.join(missing_seq)} require sequence data "
                        f"(pass --fasta)")
            changes = run_apply(doc, ctx, names)
        except KeyError as e:
            ap.error(str(e).strip('"'))
        out_text = write(doc)
        if args.out:
            with open(args.out, "w", encoding="ascii") as fh:
                fh.write(out_text)
        else:
            sys.stdout.write(out_text)
        report = render_changes(changes)
        if args.report:
            with open(args.report, "w", encoding="ascii") as fh:
                fh.write(report)
        else:
            sys.stderr.write(report)
        return 0

    ap.error("nothing to do: pass --list, --detect, or --apply")
    return 2
