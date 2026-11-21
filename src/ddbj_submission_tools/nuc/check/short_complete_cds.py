"""List complete CDS features shorter than a given amino-acid length.

An MSS annotation file (.ann) is a 5-column TSV:
    Entry / Feature / Location / Qualifier key / Qualifier value

A CDS is treated as "complete" when its location contains neither "<" nor ">"
(no partial-end markers). The coding length is the sum of all spans in the
location; the protein length is (coding length / 3) - 1 when the CDS ends with
a stop codon, so by default the stop codon is assumed to be included.

Usage:
    python -m ddbj_submission_tools.nuc.check.short_complete_cds Dioth.ann [--min-aa 30] [--no-stop]
                                       [--out short_cds.tsv]
"""

import argparse
import re
import sys

# Matches a single span such as "123..456", "123", "complement(123..456)"
SPAN_RE = re.compile(r"(\d+)\.\.(\d+)|(\d+)")


def parse_location(loc):
    """Return (length, is_complete) for an INSDC location string.

    is_complete is False when a partial marker ("<" or ">") is present.
    """
    complete = "<" not in loc and ">" not in loc
    # Strip operators and keep only the numeric spans.
    body = re.sub(r"[a-zA-Z_]+\(|\)", "", loc)
    total = 0
    for m in SPAN_RE.finditer(body):
        if m.group(1) is not None:
            start, end = int(m.group(1)), int(m.group(2))
            total += abs(end - start) + 1
        else:
            total += 1
    return total, complete


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("ann", help="MSS annotation file (.ann)")
    ap.add_argument("--min-aa", type=int, default=30,
                    help="report CDS shorter than this amino-acid length")
    ap.add_argument("--no-stop", action="store_true",
                    help="location does NOT include the stop codon")
    ap.add_argument("--out", default="-", help="output TSV ('-' for stdout)")
    args = ap.parse_args(argv)

    out = sys.stdout if args.out == "-" else open(args.out, "w")

    entry = ""
    feature = None  # dict for the CDS currently being read
    results = []
    n_cds = 0

    def flush(feat):
        """Evaluate a finished CDS feature and keep it if it is short."""
        nonlocal n_cds
        if feat is None:
            return
        n_cds += 1
        nt, complete = feat["len"], feat["complete"]
        if not complete:
            return
        aa = nt // 3 if args.no_stop else (nt // 3) - 1
        if aa < args.min_aa:
            results.append((feat["entry"], feat["loc"], nt, aa,
                            feat.get("product", ""),
                            feat.get("note", "")))

    with open(args.ann, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            cols = line.split("\t")
            cols += [""] * (5 - len(cols)) if len(cols) < 5 else []
            ent, feat_name, loc, key, val = cols[0], cols[1], cols[2], cols[3], cols[4]

            if ent:
                entry = ent
            if feat_name:  # a new feature starts here
                flush(feature)
                feature = None
                if feat_name == "CDS":
                    nt, complete = parse_location(loc)
                    feature = {"entry": entry, "loc": loc,
                               "len": nt, "complete": complete}
            elif feature is not None and key in ("product", "note"):
                feature.setdefault(key, val)
        flush(feature)

    print("#entry\tlocation\tnt_length\taa_length\tproduct\tnote", file=out)
    for row in results:
        print("\t".join(str(x) for x in row), file=out)
    print(f"# total CDS: {n_cds}; complete and < {args.min_aa} aa: "
          f"{len(results)}", file=sys.stderr)

    if out is not sys.stdout:
        out.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
