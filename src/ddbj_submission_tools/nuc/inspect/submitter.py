"""Print SUBMITTER values (ab_name or email) of the COMMON entry of an .ann file.

ab_name may appear several times; values are printed one per line in file
order (the author order on the flat file). REFERENCE values are not included.

Usage:
    python -m ddbj_submission_tools.nuc.inspect.submitter FILE.ann --field ab_name|email
"""
import argparse
import sys
from pathlib import Path

from ddbj_submission_tools.nuc import mssann


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ann", type=Path)
    ap.add_argument("--field", choices=["ab_name", "email"], required=True)
    args = ap.parse_args(argv)
    if not args.ann.is_file():
        print(f"Error: file not found: {args.ann}", file=sys.stderr)
        return 1
    try:
        common = mssann.read_common(args.ann)    # stops before the sequence entries
    except mssann.AnnParseError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    values = common.get_values("SUBMITTER", args.field) if common else []
    if not values:
        print(f"Error: COMMON/SUBMITTER/{args.field} not found in {args.ann}", file=sys.stderr)
        return 1
    for v in values:
        print(v)
    return 0


if __name__ == "__main__":
    sys.exit(main())
