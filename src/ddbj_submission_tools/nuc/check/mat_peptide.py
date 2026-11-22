"""Check mat_peptide features of an entry against the CDS they belong to.

Viral polyprotein submissions (coronaviruses, flaviviruses, bunyaviruses)
carry many mat_peptide features that are easy to get wrong by one base.
This script checks the points the DDBJ validator does not report together:

  error    the location has intervals that overlap or run backwards
               (DDBJ ANN2020). A join that repeats one base is correct
               when the parent CDS repeats the same base and carries
               /ribosomal_slippage: the mat_peptide follows the CDS, and
               /ribosomal_slippage itself cannot be given on a
               mat_peptide. That case is reported as info, because the
               ANN2020 of the DDBJ validator is a false positive there
  error    the mat_peptide is not inside any CDS of the entry
  warning  the length is not a multiple of 3
  warning  the mat_peptide is not in frame with its CDS
  warning  consecutive mat_peptides of one CDS leave a gap or overlap.
               A -1 frameshift makes two products share the slippage site
               (nsp11 and nsp12 in coronaviruses), so an overlap there is
               expected and can be ignored
  info     the mat_peptides do not cover the whole CDS (stop codon
               excluded)

A mat_peptide is assigned to the CDS with the smallest span that contains
it and has the same strand.

Usage:
    python -m ddbj_submission_tools.nuc.check.mat_peptide FILE.ann ENTRY_NAME [--fasta FILE.fasta]
                                 [--quiet]

Exit code: 0 if nothing was reported above info level, 1 if a warning or
an error was reported, 2 on a file or parse error.
"""

import argparse
import sys
from pathlib import Path

from ddbj_submission_tools.nuc import mssann
from ddbj_submission_tools.nuc.mssann import fasta


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Check mat_peptide features against their CDS.")
    parser.add_argument("ann", type=Path, help="annotation file (.ann)")
    parser.add_argument("entry", help="entry name")
    parser.add_argument("--fasta", type=Path, default=None,
                        help="sequence file (default: same base name)")
    parser.add_argument("--no-fasta-check", action="store_true",
                        help="do not read the .fasta file")
    parser.add_argument("--quiet", action="store_true",
                        help="do not print the info lines")
    return parser.parse_args(argv)


def get_seq_length(args):
    """Return (ok, seq_length); seq_length is None when the check is off."""
    if args.no_fasta_check:
        return True, None
    path = args.fasta or args.ann.with_suffix(".fasta")
    if not path.is_file():
        print(f"Error: sequence file not found: {path} "
              "(use --fasta or --no-fasta-check)", file=sys.stderr)
        return False, None
    try:
        rec = fasta.find(path, args.entry)
    except fasta.FastaParseError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return False, None
    if rec is None:
        print(f"Error: entry '{args.entry}' not found in {path}",
              file=sys.stderr)
        return False, None
    return True, rec.length


def load(feature_key, entry, seq_length, circular, problems):
    """Return [(feature, location)] for one feature key, skipping bad ones."""
    out = []
    for feat in entry.get_features(feature_key):
        if feat.location is None:
            problems.append(("error", feat.line_no,
                             f"{feat.key} has no location"))
            continue
        try:
            loc = mssann.parse_location(feat.location, seq_length=seq_length,
                                        circular=circular)
        except mssann.LocationError as exc:
            problems.append(("error", feat.line_no, str(exc)))
            continue
        out.append((feat, loc))
    return out


def label(feat):
    """Short name of a feature for the report."""
    for key in ("product", "gene", "locus_tag", "note"):
        value = feat.get(key)
        if value:
            return value
    return feat.key


def slippage_bases(cds_feat, cds_loc):
    """Return the bases a CDS repeats because of a ribosomal slippage."""
    if not cds_feat.has("ribosomal_slippage"):
        return set()
    repeated = set()
    for prev, cur in zip(cds_loc.spans, cds_loc.spans[1:]):
        if cur.start <= prev.end:
            repeated.update(range(cur.start, prev.end + 1))
    return repeated


def check_intervals(feat, loc, parent, problems):
    """Intervals must be disjoint, except at a ribosomal slippage."""
    allowed = slippage_bases(*parent) if parent is not None else set()
    for prev, cur in zip(loc.spans, loc.spans[1:]):
        if cur.start > prev.end:
            continue
        overlap = set(range(cur.start, prev.end + 1))
        if overlap and overlap <= allowed:
            shown = ", ".join(str(b) for b in sorted(overlap))
            problems.append(
                ("info", feat.line_no,
                 f"{label(feat)}: base {shown} is repeated, as in the CDS "
                 f"at line {parent[0].line_no} with /ribosomal_slippage; "
                 "correct, and the ANN2020 of the validator is a false "
                 "positive"))
        else:
            problems.append(
                ("error", feat.line_no,
                 f"{label(feat)}: overlapping location intervals "
                 f"({prev.start}..{prev.end} and {cur.start}..{cur.end}); "
                 "a repeated base is correct only where the parent CDS "
                 "repeats the same base with /ribosomal_slippage "
                 "(ANN2020)"))
        return


def assign_cds(loc, cds_list):
    """Return the smallest CDS location that contains loc, or None."""
    best = None
    for feat, cds_loc in cds_list:
        if cds_loc.strand != loc.strand and cds_loc.strand != 0:
            continue
        if cds_loc.start <= loc.start and loc.end <= cds_loc.end:
            if best is None or cds_loc.length < best[1].length:
                best = (feat, cds_loc)
    return best


def codon_start(feat):
    """Return the codon_start of a CDS as an int (1 when absent or bad)."""
    try:
        value = int(feat.get("codon_start", "1"))
    except ValueError:
        return 1
    return value if value in (1, 2, 3) else 1


def cds_offset(cds_loc, cds_feat, position):
    """Return the 0-based offset of a genomic position inside the CDS.

    Spans are walked in product order, so the value already accounts for
    introns and for the strand. None if the position is in no span.
    """
    offset = 0
    for span in cds_loc.product_spans():
        if span.start <= position <= span.end:
            inside = (position - span.start if cds_loc.strand >= 0
                      else span.end - position)
            return offset + inside - (codon_start(cds_feat) - 1)
        offset += span.length
    return None


def check_frame(feat, loc, cds_feat, cds_loc, problems):
    """The 5' end of the mat_peptide must sit on a codon boundary."""
    five_prime = loc.start if loc.strand >= 0 else loc.end
    offset = cds_offset(cds_loc, cds_feat, five_prime)
    if offset is None:
        problems.append(("warning", feat.line_no,
                         f"{label(feat)}: position {five_prime} is not "
                         f"inside the CDS at line {cds_feat.line_no}"))
    elif offset % 3 != 0:
        problems.append(("warning", feat.line_no,
                         f"{label(feat)}: not in frame with the CDS at line "
                         f"{cds_feat.line_no} (offset {offset} nt, "
                         f"{offset % 3} nt off a codon boundary)"))


def check_series(cds_feat, cds_loc, members, problems, quiet):
    """Gaps, overlaps and total coverage of one CDS."""
    members.sort(key=lambda item: item[1].start)
    for (prev_f, prev_l), (cur_f, cur_l) in zip(members, members[1:]):
        if cur_l.start > prev_l.end + 1:
            problems.append(("warning", cur_f.line_no,
                             f"{label(cur_f)}: gap of "
                             f"{cur_l.start - prev_l.end - 1} nt after "
                             f"{label(prev_f)} ({prev_l.end})"))
        elif cur_l.start <= prev_l.end:
            problems.append(("warning", cur_f.line_no,
                             f"{label(cur_f)}: overlaps {label(prev_f)} by "
                             f"{prev_l.end - cur_l.start + 1} nt"))
    if quiet:
        return
    covered = sum(loc.length for _, loc in members)
    coding = cds_loc.length - 3          # the stop codon has no peptide
    if covered != coding:
        problems.append(("info", cds_feat.line_no,
                         f"{label(cds_feat)}: mat_peptides cover {covered} nt "
                         f"of {coding} nt ({coding - covered} nt left over; "
                         "a repeated base of a frameshift counts twice)"))


def main(argv=None):
    args = parse_args(argv)
    if not args.ann.is_file():
        print(f"Error: file not found: {args.ann}", file=sys.stderr)
        return 2

    ok, seq_length = get_seq_length(args)
    if not ok:
        return 2

    try:
        common, entry = mssann.find_entry(args.ann, args.entry)
    except mssann.AnnParseError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    if entry is None:
        print(f"Error: entry '{args.entry}' not found in {args.ann}",
              file=sys.stderr)
        return 2

    circular = mssann.is_circular(entry, common)
    problems = []
    cds_list = load("CDS", entry, seq_length, circular, problems)
    peptides = load("mat_peptide", entry, seq_length, circular, problems)

    by_cds = {}
    for feat, loc in peptides:
        if loc.length % 3 != 0:
            problems.append(("warning", feat.line_no,
                             f"{label(feat)}: length {loc.length} nt is not "
                             "a multiple of 3"))
        parent = assign_cds(loc, cds_list)
        check_intervals(feat, loc, parent, problems)
        if parent is None:
            problems.append(("error", feat.line_no,
                             f"{label(feat)}: not inside any CDS of entry "
                             f"'{entry.name}'"))
            continue
        check_frame(feat, loc, parent[0], parent[1], problems)
        by_cds.setdefault(parent[0].line_no, (parent, []))[1].append(
            (feat, loc))

    for (cds_feat, cds_loc), members in by_cds.values():
        check_series(cds_feat, cds_loc, members, problems, args.quiet)

    n_bad = 0
    for level, line_no, text in sorted(problems, key=lambda p: p[1]):
        if level == "info" and args.quiet:
            continue
        print(f"{level.capitalize()}: line {line_no}: {text}",
              file=sys.stderr)
        if level != "info":
            n_bad += 1

    print(f"# entry '{entry.name}': {len(cds_list)} CDS, {len(peptides)} "
          f"mat_peptide, {n_bad} problem(s)", file=sys.stderr)
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
