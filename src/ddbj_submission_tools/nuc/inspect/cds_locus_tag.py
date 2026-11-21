"""List CDS features of an entry with their location and locus_tag.

Output is tab-separated with a header line:
    feature  location  locus_tag
With --details, parsed location columns are added:
    strand  spans  length  partial  line

Each location is checked with mssann.parse_location; problems are
reported on stderr (the row is still printed). A CDS without locus_tag
is printed with "-" in the locus_tag column.

Usage:
    python -m ddbj_submission_tools.nuc.inspect.cds_locus_tag FILE.ann ENTRY_NAME [--feature CDS]
                                 [--details] [--fasta FILE.fasta]
"""

import argparse
import sys
from pathlib import Path

from ddbj_submission_tools.nuc import mssann
from ddbj_submission_tools.nuc.mssann import fasta

MISSING = "-"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="List CDS location and locus_tag for an entry.")
    parser.add_argument("ann", type=Path, help="annotation file (.ann)")
    parser.add_argument("entry", help="entry name (e.g. chr01)")
    parser.add_argument("--feature", default="CDS",
                        help="feature key to list (default: CDS)")
    parser.add_argument("--details", action="store_true",
                        help="add parsed location columns")
    parser.add_argument("--fasta", type=Path, default=None,
                        help="sequence file (default: same base name as .ann)")
    parser.add_argument("--no-fasta-check", action="store_true",
                        help="do not read the .fasta file "
                             "(entry check and length check are skipped)")
    return parser.parse_args(argv)


def get_seq_length(args):
    """Return (ok, seq_length) after checking the entry in the .fasta file."""
    if args.no_fasta_check:
        return True, None
    fasta_path = args.fasta or args.ann.with_suffix(".fasta")
    if not fasta_path.is_file():
        print(f"Error: sequence file not found: {fasta_path} "
              "(use --fasta or --no-fasta-check)", file=sys.stderr)
        return False, None
    try:
        rec = fasta.find(fasta_path, args.entry)
    except fasta.FastaParseError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return False, None
    if rec is None:
        print(f"Error: entry '{args.entry}' not found in {fasta_path}",
              file=sys.stderr)
        return False, None
    return True, rec.length


def main(argv=None):
    args = parse_args(argv)

    if not args.ann.is_file():
        print(f"Error: file not found: {args.ann}", file=sys.stderr)
        return 1

    ok, seq_length = get_seq_length(args)
    if not ok:
        return 1

    try:
        common, entry = mssann.find_entry(args.ann, args.entry)
    except mssann.AnnParseError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    if entry is None:
        print(f"Error: entry '{args.entry}' not found in {args.ann}",
              file=sys.stderr)
        return 1

    circular = mssann.is_circular(entry, common)
    features = entry.get_features(args.feature)

    header = ["feature", "location", "locus_tag"]
    if args.details:
        header += ["strand", "spans", "length", "partial", "line"]
    print("\t".join(header))

    n_warn = 0
    for feat in features:
        tags = feat.get_all("locus_tag")
        if not tags:
            print(f"Warning: line {feat.line_no}: {feat.key} has no locus_tag",
                  file=sys.stderr)
            n_warn += 1
        elif len(tags) > 1:
            print(f"Warning: line {feat.line_no}: {feat.key} has "
                  f"{len(tags)} locus_tag values", file=sys.stderr)
            n_warn += 1
        locus_tag = tags[0] if tags else MISSING

        loc = None
        if feat.location is None:
            print(f"Warning: line {feat.line_no}: {feat.key} has no location",
                  file=sys.stderr)
            n_warn += 1
        else:
            try:
                loc = mssann.parse_location(feat.location,
                                            seq_length=seq_length,
                                            circular=circular)
            except mssann.LocationError as exc:
                print(f"Warning: line {feat.line_no}: {exc}", file=sys.stderr)
                n_warn += 1

        row = [feat.key, feat.location or MISSING, locus_tag]
        if args.details:
            if loc is None:
                row += [MISSING] * 4
            else:
                strand = {1: "+", -1: "-", 0: "mixed"}[loc.strand]
                row += [strand, str(len(loc.spans)), str(loc.length),
                        "yes" if loc.is_partial else "no"]
            row.append(str(feat.line_no))
        print("\t".join(row))

    print(f"# {len(features)} {args.feature} feature(s) in entry "
          f"'{args.entry}', {n_warn} warning(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
