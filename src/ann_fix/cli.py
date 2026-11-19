"""ann_fix: list ddbj-validator findings per code and apply chosen fixes on the .ann."""
from __future__ import annotations

import argparse
import json
import os
import sys

from .annfile import AnnFile
from .layout import merge_previous
from .short import GENOME_MIN_LENGTH, drop_short, short_report
from .operations import ARTIFICIAL_LOCATION_VALUES, DEFAULTS, IGNORABLE, KEEP_VALUES, OPERATIONS, apply, op_drop
from .report import Finding, group_by_code, parse_report


def _default_label(code: str) -> str:
    if code in DEFAULTS:
        return DEFAULTS[code][0]
    if code in IGNORABLE:
        return f"ignore ({IGNORABLE[code]})"
    return "none"


def _examples(ann: AnnFile, findings: list[Finding], n: int = 5) -> list[dict]:
    out = []
    for f in findings[:n]:
        feat = ann.feature_at_line(f.line) if f.line is not None else None
        out.append({"line": f.line, "key": feat.key if feat else f.feature,
                    "location": feat.location if feat else None,
                    "locus_tag": feat.get(ann.rows, "locus_tag") if feat else None})
    return out


def _parse_apply(spec: str, present: list[str]) -> dict[str, tuple[str, dict]]:
    if spec == "defaults":
        return {c: DEFAULTS[c] for c in present if c in DEFAULTS}
    choice: dict[str, tuple[str, dict]] = {}
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        code, eq, rest = item.partition("=")
        if not eq or not rest:
            raise ValueError(f"malformed --apply item {item!r} (expected CODE=OP[:ARG])")
        op, _, arg = rest.partition(":")
        if op not in OPERATIONS:
            raise ValueError(f"unknown operation {op!r}; choose from {sorted(OPERATIONS)}")
        kwargs = dict(DEFAULTS.get(code, ("", {}))[1]) if DEFAULTS.get(code, ("", {}))[0] == op else {}
        if arg:
            if op.startswith("overlap"):
                keep = {"trna": "other"}.get(arg.lower(), arg.lower())
                if keep not in KEEP_VALUES:
                    raise ValueError(f"keep must be one of cds/other (got {arg!r})")
                kwargs["keep"] = keep
            elif op == "artificial-location":
                if arg not in ARTIFICIAL_LOCATION_VALUES:
                    raise ValueError(f"artificial_location value must be one of {sorted(ARTIFICIAL_LOCATION_VALUES)}")
                kwargs["value"] = arg
            else:
                raise ValueError(f"{op} takes no :ARG")
        choice[code] = (op, kwargs)
    return choice


AT_OPS = {"merge-previous", "drop"}


def _run_short(ann: AnnFile, args) -> int:
    if not args.fasta:
        sys.stderr.write("error: --fasta is required with --short-seqs / --drop-short\n")
        return 2
    if args.drop_short and not (args.out and args.fasta_out):
        sys.stderr.write("error: --drop-short needs --out and --fasta-out\n")
        return 2
    r = short_report(ann, args.fasta, args.min_length)
    if args.json:
        print(json.dumps(r, ensure_ascii=False, indent=2))
    else:
        pct = 100 * r["short_bases"] / r["total_bases"] if r["total_bases"] else 0
        lengths = [s["length"] for s in r["short"]]
        print(f"{r['count']} of {r['sequences']} sequences are shorter than {r['min_length']} bp"
              + (f" ({min(lengths)}-{max(lengths)} bp, {r['short_bases']} bp, {pct:.2f}% of {r['total_bases']} bp)"
                 if lengths else ""))
        if r["lost_features"]:
            print("removing them loses: " + ", ".join(f"{k} {n}" for k, n in r["lost_features"].items()))
        if r["lost_trna_products"]:
            print("every copy on a short sequence: " + ", ".join(r["lost_trna_products"]))
    if args.drop_short:
        text = drop_short(ann, args.fasta, args.min_length)
        ann.write(args.out)
        with open(args.fasta_out, "w", newline="", encoding="utf-8") as fh:
            fh.write(text)
        print(f"wrote {args.out} and {args.fasta_out}")
    return 0


def _run_at(ann: AnnFile, args) -> int:
    if args.apply or args.detect or args.report:
        sys.stderr.write("error: --at runs alone (line numbers refer to --ann as given)\n")
        return 2
    if not args.out:
        sys.stderr.write("error: --out is required with --at\n")
        return 2
    todo = []
    for item in args.at:
        line, eq, op = item.partition("=")
        if not eq or not line.strip().isdigit() or op not in AT_OPS:
            sys.stderr.write(f"error: --at {item!r}: expected LINE=OP with OP in {sorted(AT_OPS)}\n")
            return 2
        n = int(line)
        if not 1 <= n <= len(ann.rows):
            sys.stderr.write(f"error: --at {item!r}: line out of range; nothing written\n")
            return 1
        todo.append((ann.rows[n - 1], n, op))
    log = []
    try:
        for row, n, op in sorted(todo, key=lambda t: -t[1]):     # bottom-up: earlier lines stay valid
            cur = next(i for i, r in enumerate(ann.rows) if r is row) + 1
            if op == "merge-previous":
                log.append(merge_previous(ann, cur))
            else:
                feat = ann.feature_at_line(cur)
                if feat is None:
                    raise ValueError(f"line {n} is not the first row of a feature")
                log.append(f"line {n}: {op_drop(ann, feat)}")
    except ValueError as exc:
        sys.stderr.write(f"error: {exc}; nothing written\n")
        return 1
    for line in log:
        print(line)
    ann.write(args.out)
    print(f"wrote {args.out}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="ann_fix",
                                 description="Resolve ddbj-validator findings directly on a gff2mss .ann")
    ap.add_argument("--ann", required=True)
    ap.add_argument("--report", help="validation_report.json (ddbj-validator --json) or validation_report_details.txt")
    ap.add_argument("--detect", action="store_true", help="list findings per code (default when --apply is absent)")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--apply", metavar="SPEC",
                    help="'defaults' or CODE=OP[:ARG],... (OP: misc | artificial-location | drop | overlap-drop | overlap-misc | collapse-repeats)")
    ap.add_argument("--at", metavar="LINE=OP", action="append",
                    help="apply OP to the feature starting at LINE, without a report (OP: merge-previous | drop); "
                         "for decisions ann_fix leaves to a person, e.g. the review lines of collapse-repeats")
    ap.add_argument("--out", help="output .ann (required with --apply / --at / --drop-short; may equal --ann)")
    ap.add_argument("--fasta", help="FASTA of --ann (for --short-seqs / --drop-short)")
    ap.add_argument("--short-seqs", action="store_true",
                    help="report sequences shorter than --min-length (e.g. ANN5045) and the annotation removing them loses")
    ap.add_argument("--drop-short", action="store_true",
                    help="remove those sequences from --ann and --fasta (the submitter's decision); needs --out and --fasta-out")
    ap.add_argument("--min-length", type=int, default=GENOME_MIN_LENGTH, help=f"default {GENOME_MIN_LENGTH} (genome sequences, ANN5045); other submissions often use 200")
    ap.add_argument("--fasta-out", help="output FASTA for --drop-short")
    args = ap.parse_args(argv)

    ann = AnnFile.read(args.ann)
    if args.at:
        return _run_at(ann, args)
    if args.short_seqs or args.drop_short:
        return _run_short(ann, args)
    if not args.report:
        sys.stderr.write("error: --report is required (except with --at)\n")
        return 2
    name = os.path.basename(args.ann)
    name = name[:-4] if name.endswith(".ann") else name   # the validator names a file set after the .ann
    with open(args.report, encoding="utf-8") as fh:
        text = fh.read()
    try:
        groups = group_by_code(parse_report(text, file_set=name))
    except ValueError as exc:
        sys.stderr.write(f"error: {exc}; validate this .ann under the same name first\n")
        return 1
    if not groups:
        others = sorted({f.file_set for f in parse_report(text) if f.file_set})
        if others:
            sys.stderr.write(f"note: no findings for {name} in this report (files: {', '.join(others)})\n")

    if args.apply is None:
        rows = [{"code": c, "level": fs[0].level, "count": len(fs), "default": _default_label(c),
                 "examples": _examples(ann, fs)} for c, fs in groups.items()]
        if args.json:
            print(json.dumps(rows, ensure_ascii=False, indent=2))
        else:
            for r in rows:
                print(f"{r['code']}  {r['level']}  n={r['count']}  default={r['default']}")
                for e in r["examples"]:
                    if e["line"] is not None:
                        print(f"  Line[{e['line']}] {e['key']} {e['location']} {e['locus_tag'] or ''}")
        return 0

    if not args.out:
        sys.stderr.write("error: --out is required with --apply\n")
        return 2
    try:
        choice = _parse_apply(args.apply, list(groups))
    except ValueError as exc:
        sys.stderr.write(f"error: {exc}\n")
        return 2
    findings = [f for fs in groups.values() for f in fs]
    try:
        log = apply(ann, findings, choice)
    except ValueError as exc:
        sys.stderr.write(f"error: {exc}; nothing written\n")
        return 1
    for line in log:
        print(line)
    ann.write(args.out)
    print(f"wrote {args.out}")
    print("re-validate: ddbj-validator ddbj <dir containing the .ann and .fasta> -o <dir>_out -f -j 1")
    return 0


if __name__ == "__main__":
    sys.exit(main())
