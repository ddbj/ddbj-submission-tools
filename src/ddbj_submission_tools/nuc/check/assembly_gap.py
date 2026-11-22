"""Check the assembly_gap features of an .ann file against the runs of N in its FASTA.

First the runs of N (at least --min-gap-length, default 10) are counted by length, and
the share of runs whose length is a multiple of 100 is used to suggest estimated_length
for the runs that have no feature:

  - mostly multiples of 100 (>= --unknown-fraction, default 0.5): the N were put in by
    a scaffolder as gaps of unknown length -> estimated_length "unknown"
  - otherwise (lengths vary): the lengths were estimated -> estimated_length "known"

The suggestion is only a starting point: gap_type and linkage_evidence come from how
the sequences were assembled, which only the submitter knows.

Then every gap / assembly_gap feature is checked:
  - its location is exactly one run of N (no join, no extra or missing bases)
  - estimated_length and gap_type are present and valid; an integer estimated_length
    equals the run length
  - linkage_evidence is present (and valid) when gap_type requires it
  - runs of N that no feature covers are listed (add them with nuc.fix.assembly_gap)
  - a FASTA sequence with runs of N but no entry in the .ann is reported

Usage:
    python -m ddbj_submission_tools.nuc.check.assembly_gap FILE.ann FILE.fasta [--min-gap-length 10]

Exit code: 0 no problem, 1 problems or runs without a feature, 2 usage or file error.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from dataclasses import dataclass, field

from ddbj_submission_tools.nuc.ann_fix.annfile import AnnFile
from ddbj_submission_tools.nuc.common.fasta import parse_fasta_sequences

GAP_TYPES = {"between scaffolds", "within scaffold", "telomere", "centromere", "short arm",
             "heterochromatin", "repeat within scaffold", "repeat between scaffolds",
             "contamination", "unknown"}
LINKAGE_EVIDENCE = {"align genus", "align xgenus", "align trnscpt", "clone contig", "map",
                    "paired-ends", "pcr", "proximity ligation", "strobe", "within clone",
                    "unspecified"}
NEEDS_LINKAGE = {"within scaffold", "repeat within scaffold", "contamination"}
GAP_KEYS = {"assembly_gap", "gap"}
_N_RUN = re.compile("[Nn]+")
_SPAN = re.compile(r"^(\d+)\.\.(\d+)$|^(\d+)$")


@dataclass
class Problem:
    entry: str
    location: str
    kind: str          # location | estimated_length | gap_type | linkage_evidence | entry
    message: str


@dataclass
class Analysis:
    lengths: list[int] = field(default_factory=list)                     # every run of N
    uncovered: dict[str, list[tuple[int, int]]] = field(default_factory=dict)
    problems: list[Problem] = field(default_factory=list)


def n_runs(seq: str, min_len: int) -> list[tuple[int, int]]:
    """1-based (start, end) of every run of N/n of at least *min_len* bases."""
    return [(m.start() + 1, m.end()) for m in _N_RUN.finditer(seq) if m.end() - m.start() >= min_len]


def suggest_estimated_length(lengths: list[int], unknown_fraction: float = 0.5) -> tuple[str, float]:
    """("unknown" | "known", share of lengths that are multiples of 100)."""
    if not lengths:
        return "known", 0.0
    frac = sum(1 for n in lengths if n % 100 == 0) / len(lengths)
    return ("unknown" if frac >= unknown_fraction else "known"), frac


def _span(location: str) -> tuple[int, int] | None:
    m = _SPAN.match(location.strip())
    if not m:
        return None
    if m.group(3):
        return int(m.group(3)), int(m.group(3))
    return int(m.group(1)), int(m.group(2))


def _check_qualifiers(ann: AnnFile, feat, problems: list[Problem], run_len: int | None) -> None:
    add = lambda kind, msg: problems.append(Problem(feat.entry, feat.location, kind, msg))
    est = feat.get(ann.rows, "estimated_length")
    if est is None:
        add("estimated_length", "no estimated_length")
    elif est not in ("known", "unknown"):
        if not est.isdigit():
            add("estimated_length", f"estimated_length '{est}' is not known, unknown or a number")
        elif run_len is not None and int(est) != run_len:
            add("estimated_length", f"estimated_length {est} but the run of N is {run_len} bases")
    if feat.key != "assembly_gap":
        return
    gap_type = feat.get(ann.rows, "gap_type")
    evidence = [v for k, v in feat.qualifiers(ann.rows) if k == "linkage_evidence"]
    if gap_type is None:
        add("gap_type", "no gap_type")
    elif gap_type not in GAP_TYPES:
        add("gap_type", f"gap_type '{gap_type}' is not a valid value")
    if gap_type in NEEDS_LINKAGE and not evidence:
        add("linkage_evidence", f"gap_type '{gap_type}' needs linkage_evidence")
    for v in evidence:
        if v not in LINKAGE_EVIDENCE:
            add("linkage_evidence", f"linkage_evidence '{v}' is not a valid value")


def analyze(ann: AnnFile, seqs: dict[str, str], min_len: int = 10) -> Analysis:
    result = Analysis()
    feats_by_entry: dict[str, list] = {}
    for f in ann.features():
        feats_by_entry.setdefault(f.entry, []).append(f)

    for name, s in seqs.items():
        runs = n_runs(s, min_len)
        result.lengths += [e - b + 1 for b, e in runs]
        if name not in feats_by_entry:
            if runs:
                result.problems.append(Problem(name, "", "entry", f"{len(runs)} run(s) of N but no entry in the .ann"))
            continue
        run_set = set(runs)
        covered: list[tuple[int, int]] = []
        for f in feats_by_entry[name]:
            if f.key not in GAP_KEYS:
                continue
            span = _span(f.location)
            if span is None or span not in run_set:
                result.problems.append(Problem(name, f.location, "location",
                                               f"{f.key} {f.location} is not exactly one run of N"))
            if span is not None:
                covered.append(span)
            _check_qualifiers(ann, f, result.problems, span[1] - span[0] + 1 if span in run_set else None)
        missing = [r for r in runs if not any(b <= r[1] and r[0] <= e for b, e in covered)]
        if missing:
            result.uncovered[name] = missing
    return result


def report(result: Analysis, unknown_fraction: float, out=None) -> None:
    out = out or sys.stdout
    if not result.lengths:
        print("# no runs of N", file=out)
        for p in result.problems:
            print(f"ERROR\t{p.entry}\t{p.location}\t{p.kind}\t{p.message}", file=out)
        return
    counts = Counter(result.lengths)
    est, frac = suggest_estimated_length(result.lengths, unknown_fraction)
    n_missing = sum(len(v) for v in result.uncovered.values())
    print(f"# runs of N: {len(result.lengths)} ({len(counts)} distinct lengths); "
          f"multiples of 100: {frac:.0%}", file=out)
    print("length\tcount", file=out)
    for length, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:20]:
        print(f"{length}\t{n}", file=out)
    if len(counts) > 20:
        print(f"... {len(counts) - 20} other lengths", file=out)
    if est == "unknown":
        print(f"# suggestion: estimated_length unknown ({frac:.0%} of the runs are multiples of 100: "
              "gaps of unknown length put in while scaffolding)", file=out)
    else:
        print(f"# suggestion: estimated_length known (only {frac:.0%} of the runs are multiples of 100: "
              "the lengths vary, so they were estimated)", file=out)
    print(f"# runs without a gap feature: {n_missing}", file=out)
    for name, runs in result.uncovered.items():
        for b, e in runs:
            print(f"MISSING\t{name}\t{b}..{e}\t{e - b + 1}", file=out)
    for p in result.problems:
        print(f"ERROR\t{p.entry}\t{p.location}\t{p.kind}\t{p.message}", file=out)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ann", help="annotation file (.ann)")
    ap.add_argument("fasta", help="sequence file (FASTA, may end entries with //)")
    ap.add_argument("--min-gap-length", type=int, default=10, help="shortest run of N treated as a gap (default 10)")
    ap.add_argument("--unknown-fraction", type=float, default=0.5,
                    help="share of multiples of 100 from which 'unknown' is suggested (default 0.5)")
    return ap.parse_args(argv)


def load(ann_path: str, fasta_path: str) -> tuple[AnnFile, dict[str, str]]:
    return AnnFile.read(ann_path), parse_fasta_sequences(fasta_path)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        ann, seqs = load(args.ann, args.fasta)
    except (OSError, UnicodeDecodeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2
    result = analyze(ann, seqs, args.min_gap_length)
    report(result, args.unknown_fraction)
    return 1 if result.problems or result.uncovered else 0


if __name__ == "__main__":
    sys.exit(main())
