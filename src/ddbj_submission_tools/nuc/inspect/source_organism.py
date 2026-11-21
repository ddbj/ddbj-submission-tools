"""Print the organism of the source feature for a given entry name.

The entry name is the value in the Entry column of the .ann file (other
than COMMON) and corresponds to the ">" line of the .fasta file.

If the entry has no source feature of its own, the COMMON source
(location "1..E") is used, as defined in the spec.

Usage:
    python -m ddbj_submission_tools.nuc.inspect.source_organism FILE.ann ENTRY_NAME [--fasta FILE.fasta]

By default the .fasta file with the same base name as FILE.ann is used to
check that ENTRY_NAME exists in the sequence file.
"""

import argparse
import sys
from pathlib import Path

from ddbj_submission_tools.nuc import mssann
from ddbj_submission_tools.nuc.mssann import fasta


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Print source/organism for an entry of an MSS .ann file.")
    parser.add_argument("ann", type=Path, help="annotation file (.ann)")
    parser.add_argument("entry", help="entry name (e.g. sequence1)")
    parser.add_argument("--fasta", type=Path, default=None,
                        help="sequence file (default: same base name as .ann)")
    parser.add_argument("--no-fasta-check", action="store_true",
                        help="skip checking the entry name in the .fasta file")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if not args.ann.is_file():
        print(f"Error: file not found: {args.ann}", file=sys.stderr)
        return 1

    # Check that the entry name exists in the sequence file
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

    # Find the entry (and COMMON) in the annotation file
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

    organism = source.get("organism")
    if organism is None:
        print(f"Error: source of entry '{args.entry}' has no organism "
              f"(line {source.line_no})", file=sys.stderr)
        return 1

    if entry is None or entry.get_feature("source") is None:
        print(f"Note: organism taken from COMMON source "
              f"(line {source.line_no})", file=sys.stderr)

    print(organism)
    return 0


if __name__ == "__main__":
    sys.exit(main())
