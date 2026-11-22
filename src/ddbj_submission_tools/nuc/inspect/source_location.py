"""Print the location (3rd column) of the source feature for an entry name.

The location is shown as written in the .ann file, followed by its parsed
form (operator, strand, spans with partial markers, total length).

If the entry has no source of its own, the COMMON source is used and the
meta position "E" is replaced with the sequence length from the .fasta file.

Usage:
    python -m ddbj_submission_tools.nuc.inspect.source_location FILE.ann ENTRY_NAME [--fasta FILE.fasta]
    python -m ddbj_submission_tools.nuc.inspect.source_location FILE.ann ENTRY_NAME --raw
"""

import argparse
import sys
from pathlib import Path

from ddbj_submission_tools.nuc import mssann
from ddbj_submission_tools.nuc.mssann import fasta


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Print source location for an entry of an MSS .ann file.")
    parser.add_argument("ann", type=Path, help="annotation file (.ann)")
    parser.add_argument("entry", help="entry name (e.g. sequence1)")
    parser.add_argument("--fasta", type=Path, default=None,
                        help="sequence file (default: same base name as .ann)")
    parser.add_argument("--no-fasta-check", action="store_true",
                        help="do not read the .fasta file "
                             "(entry check, 'E' and length check are skipped)")
    parser.add_argument("--raw", action="store_true",
                        help="print only the location string as written")
    return parser.parse_args(argv)


def yes_no(flag):
    return "yes" if flag else "no"


def print_details(entry_name, source, loc, seq_length, from_common, circular):
    """Print the parsed location in a human-readable form."""
    origin = "COMMON" if from_common else entry_name
    print(f"entry      : {entry_name}")
    print(f"feature    : source (defined in {origin}, line {source.line_no})")
    print(f"location   : {source.location}")
    if loc.raw != source.location:
        print(f"resolved   : {loc.raw}")
    print(f"operator   : {loc.operator or '-'}")
    print(f"complement : {yes_no(loc.complement)}")
    strand = {1: "+", -1: "-", 0: "mixed"}[loc.strand]
    print(f"strand     : {strand}")
    print(f"topology   : {'circular' if circular else 'linear'}")
    print(f"range      : {loc.start}..{loc.end}")
    print(f"length     : {loc.length} bp (sum of spans)")
    print(f"spans      : {len(loc.spans)}")
    for i, span in enumerate(loc.spans, start=1):
        print(f"  [{i}] {span}  start={span.start} end={span.end} "
              f"partial_start={yes_no(span.partial_start)} "
              f"partial_end={yes_no(span.partial_end)} "
              f"complement={yes_no(span.complement)} "
              f"between={yes_no(span.between)}")
    if seq_length is not None:
        print(f"seq length : {seq_length} bp (.fasta)")


def main(argv=None):
    args = parse_args(argv)

    if not args.ann.is_file():
        print(f"Error: file not found: {args.ann}", file=sys.stderr)
        return 1

    # Look up the entry in the sequence file to get its length
    seq_length = None
    if not args.no_fasta_check:
        fasta_path = args.fasta or args.ann.with_suffix(".fasta")
        if not fasta_path.is_file():
            print(f"Error: sequence file not found: {fasta_path} "
                  "(use --fasta or --no-fasta-check)", file=sys.stderr)
            return 1
        try:
            rec = fasta.find(fasta_path, args.entry)
        except fasta.FastaParseError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
        if rec is None:
            print(f"Error: entry '{args.entry}' not found in {fasta_path}",
                  file=sys.stderr)
            return 1
        seq_length = rec.length

    try:
        common, entry = mssann.find_entry(args.ann, args.entry)
    except mssann.AnnParseError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    source = mssann.effective_source(entry, common)
    if source is None:
        print(f"Error: no source feature for entry '{args.entry}' "
              f"in {args.ann}", file=sys.stderr)
        return 1
    if source.location is None:
        print(f"Error: source has no location (line {source.line_no})",
              file=sys.stderr)
        return 1

    if args.raw:
        print(source.location)
        return 0

    # TOPOLOGY circular allows "seq_length^1" in locations
    circular = mssann.is_circular(entry, common)

    try:
        loc = mssann.parse_location(source.location, seq_length=seq_length,
                                    circular=circular)
    except mssann.LocationError as exc:
        print(f"Error: line {source.line_no}: {exc}", file=sys.stderr)
        return 1

    # Show the resolved form (e.g. "1..E" -> "1..4618897") if it differs
    if "E" in source.location:
        loc.raw = source.location.replace("E", str(seq_length))

    from_common = entry is None or entry.get_feature("source") is None
    print_details(args.entry, source, loc, seq_length, from_common, circular)

    # source must cover the whole sequence
    if seq_length is not None and (loc.start != 1 or loc.end != seq_length):
        print(f"Warning: source {loc.start}..{loc.end} does not match "
              f"sequence length {seq_length} (source must cover the whole sequence)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
