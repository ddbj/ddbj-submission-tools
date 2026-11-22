"""Add assembly_gap features for the runs of N that have no gap feature in the .ann.

Runs of N are found as in nuc.check.assembly_gap. Runs already covered by a gap or
assembly_gap feature (right or wrong) are left alone; the new features go at the end of
their entry, in the same form as mss_builder writes them.

  --estimated-length auto (default) takes the suggestion of nuc.check.assembly_gap:
  "unknown" when most runs are multiples of 100 (gaps put in while scaffolding),
  "known" when the lengths vary. Run the check first and confirm it.
  --gap-type (default "within scaffold") and --linkage-evidence describe how the
  sequences were joined; take them from the submitter, do not guess.
  linkage_evidence is required for within scaffold, repeat within scaffold and contamination.

Usage:
    python -m ddbj_submission_tools.nuc.fix.assembly_gap FILE.ann FILE.fasta \\
        --linkage-evidence paired-ends [--gap-type "within scaffold"] [--estimated-length auto] --dry-run
    ... --outdir DIR | --in-place   (--in-place keeps FILE.ann.orig)

Exit code: 0 done (or nothing to add), 1 some runs are overlapped by a wrong gap feature
(fix those by hand; see nuc.check.assembly_gap), 2 usage or file error.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from ddbj_submission_tools.nuc.check.assembly_gap import (GAP_TYPES, LINKAGE_EVIDENCE, NEEDS_LINKAGE, analyze,
                                                          load, suggest_estimated_length)
from ddbj_submission_tools.nuc.common.gap_annotator import GapAnnotator


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ann", type=Path, help="annotation file (.ann)")
    ap.add_argument("fasta", type=Path, help="sequence file (FASTA)")
    ap.add_argument("--gap-type", default="within scaffold", help="gap_type (default 'within scaffold')")
    ap.add_argument("--linkage-evidence", help="linkage_evidence (e.g. paired-ends, align genus, proximity ligation)")
    ap.add_argument("--estimated-length", choices=["auto", "known", "unknown"], default="auto")
    ap.add_argument("--min-gap-length", type=int, default=10, help="shortest run of N treated as a gap (default 10)")
    ap.add_argument("--unknown-fraction", type=float, default=0.5,
                    help="share of multiples of 100 from which auto chooses 'unknown' (default 0.5)")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="report only; write nothing")
    mode.add_argument("--outdir", type=Path, help="write the new .ann to this directory")
    mode.add_argument("--in-place", action="store_true", help="overwrite the .ann, keeping FILE.ann.orig")
    return ap.parse_args(argv)


def _usage_error(args) -> str | None:
    if args.gap_type not in GAP_TYPES:
        return f"--gap-type '{args.gap_type}' is not one of: {', '.join(sorted(GAP_TYPES))}"
    if args.linkage_evidence is not None and args.linkage_evidence not in LINKAGE_EVIDENCE:
        return f"--linkage-evidence '{args.linkage_evidence}' is not one of: {', '.join(sorted(LINKAGE_EVIDENCE))}"
    if args.gap_type in NEEDS_LINKAGE and args.linkage_evidence is None:
        return f"gap_type '{args.gap_type}' needs --linkage-evidence (ask the submitter how the sequences were joined)"
    return None


def main(argv=None) -> int:
    args = parse_args(argv)
    err = _usage_error(args)
    if err:
        print(f"Error: {err}", file=sys.stderr)
        return 2
    try:
        ann, seqs = load(str(args.ann), str(args.fasta))
    except (OSError, UnicodeDecodeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2
    result = analyze(ann, seqs, args.min_gap_length)
    est = args.estimated_length
    if est == "auto":
        est, frac = suggest_estimated_length(result.lengths, args.unknown_fraction)
        print(f"estimated_length: {est} (auto; {frac:.0%} of the runs of N are multiples of 100)", file=sys.stderr)
    annotator = GapAnnotator(args.linkage_evidence or "", min_gap_length=args.min_gap_length,
                             gap_type=args.gap_type, estimated_length=est)

    # insert from the bottom of the file up, so row indices above stay valid
    last_row = {}
    for f in ann.features():
        last_row[f.entry] = f.last
    added = 0
    for name in sorted(result.uncovered, key=lambda n: last_row[n], reverse=True):
        wanted = {f"{b}..{e}" for b, e in result.uncovered[name]}
        rows = annotator.annotate(seqs[name])
        at = last_row[name] + 1
        for i in range(0, len(rows), 3):
            loc = rows[i][2]
            if loc not in wanted:
                continue
            quals = [(r[3], r[4]) for r in rows[i:i + 3] if not (r[3] == "linkage_evidence" and not r[4])]
            feat = ann.insert_feature(at, "", "assembly_gap", loc, quals)
            at = feat.last + 1
            added += 1
            print(f"add\t{name}\t{loc}\t{est}", file=sys.stderr)

    left = [p for p in result.problems if p.kind == "location"]
    for p in left:
        print(f"left\t{p.entry}\t{p.location}\t{p.message}", file=sys.stderr)
    verb = "would add" if args.dry_run else "added"
    print(f"# {verb} {added} assembly_gap feature(s); {len(left)} wrong gap feature(s) left to fix by hand",
          file=sys.stderr)
    if added and not args.dry_run:
        if args.in_place:
            shutil.copy2(args.ann, str(args.ann) + ".orig")
            target = args.ann
        else:
            args.outdir.mkdir(parents=True, exist_ok=True)
            target = args.outdir / args.ann.name
            if target.resolve() == args.ann.resolve():
                print(f"Error: --outdir is the input directory: {args.outdir}", file=sys.stderr)
                return 2
        ann.write(target)
    return 1 if left else 0


if __name__ == "__main__":
    sys.exit(main())
