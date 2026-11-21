"""Extract the sequence of an arbitrary location from an entry.

Unlike extract_cds_seq.py / extract_trna_seq.py / extract_rrna_seq.py,
which work on feature keys, this script takes the location as an argument.
It is for looking at a region that does not correspond to a feature, or
at a feature whose boundaries you want to try changing.

The location is parsed with mssann.parse_location and cut out with
Location.extract, so the full notation is accepted:

    3784..4474                      a plain range
    complement(10350..11949)        reverse complement
    join(56696..57377,1..179)       joined spans, in product order
    complement(join(10..20,30..40)) spans reversed, then complemented
    <1..41                          partial markers are accepted and
                                    reported, and do not change the cut

With --translate the region is translated with Biopython. Use
--transl-table and --codon-start to say how; the defaults (table 1,
codon_start 1) are NOT taken from any feature, because the region need
not be one.

Output is FASTA by default. Header format:

    >ENTRY:LOCATION entry=ENTRY location=LOCATION length=N [aa=N]

With --plain the sequence is printed with coordinates in the left margin,
which is easier to read when checking a boundary by eye.

Usage:
    python -m ddbj_submission_tools.nuc.inspect.region FILE.ann ENTRY_NAME LOCATION
    python -m ddbj_submission_tools.nuc.inspect.region FILE.fasta ENTRY_NAME LOCATION
    python -m ddbj_submission_tools.nuc.inspect.region FILE.ann ENTRY LOCATION --plain
    python -m ddbj_submission_tools.nuc.inspect.region FILE.ann ENTRY LOCATION --translate \
        --transl-table 5 --codon-start 2

The first argument may be the .ann or the .fasta. With the .ann, the
sequence file is found by replacing the suffix (override with --fasta),
and TOPOLOGY circular of the entry is honoured, which allows a location
that crosses the origin (e.g. join(56696..57377,1..179)).
"""

import argparse
import sys
from pathlib import Path

from ddbj_submission_tools.nuc import mssann
from ddbj_submission_tools.nuc.mssann import fasta


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract the sequence of a location from an entry.")
    parser.add_argument("ann", type=Path,
                        help="annotation file (.ann) or sequence file (.fasta)")
    parser.add_argument("entry", help="entry name (e.g. Parapseudes_sp.)")
    parser.add_argument("location",
                        help="location (e.g. '3784..4474', "
                             "'complement(<1..41)', 'join(10..20,30..40)')")
    parser.add_argument("--fasta", type=Path, default=None,
                        help="sequence file (default: same base name as .ann)")
    parser.add_argument("--translate", action="store_true",
                        help="also output the amino acid sequence "
                             "(needs Biopython)")
    parser.add_argument("--transl-table", type=int, default=1,
                        help="genetic code for --translate (default: 1)")
    parser.add_argument("--codon-start", type=int, default=1,
                        choices=(1, 2, 3),
                        help="codon_start for --translate (default: 1)")
    parser.add_argument("--plain", action="store_true",
                        help="print with coordinates instead of FASTA")
    parser.add_argument("--width", type=int, default=60,
                        help="line width of output sequences (0: no wrap)")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="output file (default: stdout)")
    return parser.parse_args(argv)


def wrap(seq, width):
    if width <= 0:
        return seq
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def plain_lines(seq, loc, width):
    """Yield the sequence with the coordinate of each line in the margin.

    Coordinates are only shown for a single span on the plus strand; for
    anything else the position within the product is shown instead, since
    the entry coordinate does not run left to right.
    """
    simple = (len(loc.spans) == 1 and not loc.complement
              and not loc.spans[0].complement)
    origin = loc.spans[0].start if simple else 1
    label = "" if simple else "+"
    step = width if width > 0 else len(seq)
    for i in range(0, len(seq), step):
        yield f"{label}{origin + i:>7}  {seq[i:i + step]}"


def report(loc, seq_length, circular):
    """Write what was parsed to stderr, so stdout stays pure FASTA."""
    strand = {1: "+", -1: "-", 0: "mixed"}[loc.strand]
    print(f"Info: entry length {seq_length} bp, "
          f"topology {'circular' if circular else 'linear'}", file=sys.stderr)
    print(f"Info: location {loc.raw}: operator {loc.operator or '-'}, "
          f"strand {strand}, {len(loc.spans)} span(s), {loc.length} bp",
          file=sys.stderr)
    if loc.is_partial:
        marks = []
        if loc.five_prime_partial:
            marks.append("5' end")
        if loc.three_prime_partial:
            marks.append("3' end")
        print(f"Info: partial at {' and '.join(marks)}; "
              "the partial markers do not change the extracted sequence",
              file=sys.stderr)
    if loc.is_between:
        print("Info: the location contains a 'a^b' site (length 0)",
              file=sys.stderr)


def translate(seq, table, codon_start):
    try:
        from Bio.Seq import Seq
    except ImportError:
        print("Error: --translate needs Biopython", file=sys.stderr)
        return None
    trimmed = seq[codon_start - 1:]
    extra = len(trimmed) % 3
    if extra:
        print(f"Info: length is not a multiple of 3 ({len(trimmed)} bp "
              f"after codon_start {codon_start}); "
              f"trailing {extra} base(s) ignored", file=sys.stderr)
        trimmed = trimmed[:len(trimmed) - extra]
    try:
        aa = str(Seq(trimmed).translate(table=table))
    except Exception as exc:                      # invalid table number
        print(f"Error: cannot translate with table {table}: {exc}",
              file=sys.stderr)
        return None
    internal = aa[:-1].count("*") if aa.endswith("*") else aa.count("*")
    if internal:
        print(f"Warning: {internal} internal stop codon(s)", file=sys.stderr)
    if aa and not aa.endswith("*"):
        print("Info: does not end with a stop codon", file=sys.stderr)
    return aa


def main(argv=None):
    args = parse_args(argv)

    if not args.ann.is_file():
        print(f"Error: file not found: {args.ann}", file=sys.stderr)
        return 1

    # The first argument may be the .fasta itself
    is_fasta = args.ann.suffix.lower() in (".fasta", ".fa", ".fna", ".seq")
    fasta_path = args.fasta or (args.ann if is_fasta
                                else args.ann.with_suffix(".fasta"))
    if not fasta_path.is_file():
        print(f"Error: sequence file not found: {fasta_path} (use --fasta)",
              file=sys.stderr)
        return 1

    try:
        seq = fasta.read_sequence(fasta_path, args.entry)
    except fasta.FastaParseError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    if seq is None:
        print(f"Error: entry '{args.entry}' not found in {fasta_path}",
              file=sys.stderr)
        return 1

    # TOPOLOGY circular allows a location that crosses the origin
    circular = False
    if not is_fasta:
        try:
            common, entry = mssann.find_entry(args.ann, args.entry)
        except mssann.AnnParseError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
        if entry is None:
            print(f"Warning: entry '{args.entry}' not found in {args.ann}; "
                  "TOPOLOGY is treated as linear", file=sys.stderr)
        else:
            circular = mssann.is_circular(entry, common)

    try:
        loc = mssann.parse_location(args.location, seq_length=len(seq),
                                    circular=circular)
    except mssann.LocationError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    sub = loc.extract(seq)
    report(loc, len(seq), circular)

    aa = None
    if args.translate:
        aa = translate(sub, args.transl_table, args.codon_start)
        if aa is None:
            return 1

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        if args.plain:
            for line in plain_lines(sub, loc, args.width):
                out.write(line + "\n")
            if aa is not None:
                out.write("\n")
                out.write(wrap(aa, args.width) + "\n")
        else:
            header = (f">{args.entry}:{loc.raw} entry={args.entry} "
                      f"location={loc.raw} length={len(sub)}")
            out.write(header + "\n" + wrap(sub, args.width) + "\n")
            if aa is not None:
                out.write(header.replace(f"length={len(sub)}",
                                         f"length={len(sub)} aa={len(aa)}")
                          + "\n" + wrap(aa, args.width) + "\n")
    finally:
        if args.output:
            out.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
