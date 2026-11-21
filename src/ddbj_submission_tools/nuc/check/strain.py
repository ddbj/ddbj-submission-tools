"""Check strain values that start with a culture collection code.

A strain value that starts with a culture collection code from NCBI's
coll_dump.txt (column 2 contains "c") must have a space between the code
and the number ("JCM11017" -> "JCM 11017"). Codes in the separator
exception list use their own separator (NIES-1234).

Get NCBI's coll_dump.txt first and pass it with --coll-dump:
    wget https://ftp.ncbi.nih.gov/pub/taxonomy/coll_dump.txt

Exit status: 0 if all strain values are correct, 1 if any error,
2 if coll_dump.txt or the exception list cannot be used.

Usage:
    python -m ddbj_submission_tools.nuc.check.strain FILE.ann [FILE.ann ...] [--verbose] [--tsv]
        --coll-dump PATH [--exceptions PATH]
"""

import argparse
import sys
from pathlib import Path

from ddbj_submission_tools.nuc import mssann
from ddbj_submission_tools.nuc.mssann import strain

COLL_DUMP_URL = "https://ftp.ncbi.nih.gov/pub/taxonomy/coll_dump.txt"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Check culture collection numbers in strain.")
    parser.add_argument("ann", type=Path, nargs="+", help="annotation file(s)")
    parser.add_argument("--coll-dump", type=Path, default=None,
                        help=f"NCBI coll_dump.txt ({COLL_DUMP_URL})")
    parser.add_argument("--exceptions", type=Path, default=None,
                        help="separator exception list "
                             f"(default: {strain.DEFAULT_EXCEPTIONS})")
    parser.add_argument("--verbose", action="store_true",
                        help="also report correct strain values")
    parser.add_argument("--tsv", action="store_true",
                        help="print all results as TSV on stdout")
    return parser.parse_args(argv)


def iter_strains(path):
    """Yield (entry name, line number, value) for every strain qualifier."""
    for entry in mssann.parse(path):
        for feat in entry.features:
            for q in feat.qualifiers:
                if q.key == "strain":
                    yield entry.name, q.line_no, q.value


def main(argv=None):
    args = parse_args(argv)
    if args.coll_dump is None:
        print(f"Error: give NCBI's coll_dump.txt with --coll-dump (get it with: wget {COLL_DUMP_URL})",
              file=sys.stderr)
        return 2
    try:
        codes = strain.load_culture_collection_codes(args.coll_dump)
        exceptions = strain.load_separator_exceptions(args.exceptions)
    except (OSError, strain.StrainDataError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    if args.tsv:
        print("\t".join(["file", "line", "entry", "strain", "result", "code",
                         "fix"]))

    n_value = n_error = 0
    for path in args.ann:
        if not path.is_file():
            print(f"Error: file not found: {path}", file=sys.stderr)
            n_error += 1
            continue
        try:
            values = list(iter_strains(path))
        except mssann.AnnParseError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            n_error += 1
            continue
        for entry_name, line_no, value in values:
            n_value += 1
            res = strain.check_strain_value(value, codes, exceptions)
            where = f"{path}:{line_no}: {entry_name}: strain '{value}'"
            if not res.ok:
                n_error += 1
                print(f"Error: {where}: {res.message}; fix: '{res.fix}'",
                      file=sys.stderr)
            elif args.verbose:
                print(f"OK: {where}: {res.message}", file=sys.stderr)
            if args.tsv:
                print("\t".join([str(path), str(line_no), entry_name, value,
                                 "OK" if res.ok else "ERROR", res.code or "",
                                 res.fix or ""]))

    print(f"# {n_value} strain value(s) checked in {len(args.ann)} file(s), "
          f"{n_error} error(s); {len(codes)} culture collection codes",
          file=sys.stderr)
    return 1 if n_error else 0


if __name__ == "__main__":
    sys.exit(main())
