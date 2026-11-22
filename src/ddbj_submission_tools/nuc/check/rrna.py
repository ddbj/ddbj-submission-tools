"""Extract rRNA sequences of an entry as FASTA and check their naming.

For each rRNA feature, the location is parsed and the sequence is cut out
of the .fasta entry (join and complement handled; the output is on the
rRNA's own strand, 5' -> 3').

Output header:
    >ID entry=ENTRY feature=rRNA location=LOCATION length=N
        [gene=...] [product=...] [pseudogene=VALUE | pseudo]
ID is the locus_tag, else the gene name (made unique with _2, _3, ...).

Checks (stderr):
  Warning: no product; product not in a standard form ("16S ribosomal RNA",
           "16S rRNA", "small/large subunit ribosomal RNA"); size in the
           gene name (rrn16, rrs, ...) and the product disagree;
           location errors; pseudo / pseudogene problems.
  Info:    partial location ('<' / '>'); missing locus_tag and quoted
           locations (once per file).

Usage:
    python -m ddbj_submission_tools.nuc.check.rrna FILE.ann ENTRY_NAME [--fasta FILE.fasta]
        [--rna] [--width 60] [-o OUT.fasta]
"""

import argparse
import sys
from pathlib import Path

from ddbj_submission_tools.nuc import mssann
from ddbj_submission_tools.nuc.mssann import fasta, rrna
from ddbj_submission_tools.nuc.mssann.seqfeature import SeqIdMaker, iter_feature_seqs, partial_note

FEATURE = "rRNA"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract rRNA sequences of an entry as FASTA.")
    parser.add_argument("ann", type=Path, help="annotation file (.ann)")
    parser.add_argument("entry", help="entry name (e.g. cp3Moroko)")
    parser.add_argument("--fasta", type=Path, default=None,
                        help="sequence file (default: same base name as .ann)")
    parser.add_argument("--rna", action="store_true",
                        help="output RNA letters (T -> U)")
    parser.add_argument("--width", type=int, default=60,
                        help="line width of output sequences (0: no wrap)")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="output file (default: stdout)")
    return parser.parse_args(argv)


def wrap(seq, width):
    if width <= 0:
        return seq
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def check_rrna(feat, warn):
    """Check product form and gene / product size agreement."""
    gene = feat.get("gene")
    product = feat.get("product")
    if product is None:
        warn("no product")
        return
    product_class = rrna.parse_rrna_product(product)
    if product_class is None:
        warn(f"product '{product}' is not in a standard form "
             "(e.g. '16S ribosomal RNA', 'small subunit ribosomal RNA')")
        return
    if gene is not None:
        gene_size = rrna.rrna_gene_size(gene)
        if gene_size is not None and product_class not in ("SSU", "LSU") \
                and gene_size != product_class:
            warn(f"gene '{gene}' ({gene_size}) and product '{product}' "
                 "disagree")


def main(argv=None):
    args = parse_args(argv)

    if not args.ann.is_file():
        print(f"Error: file not found: {args.ann}", file=sys.stderr)
        return 1
    fasta_path = args.fasta or args.ann.with_suffix(".fasta")
    if not fasta_path.is_file():
        print(f"Error: sequence file not found: {fasta_path} (use --fasta)",
              file=sys.stderr)
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

    try:
        seq = fasta.read_sequence(fasta_path, args.entry)
    except fasta.FastaParseError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    if seq is None:
        print(f"Error: entry '{args.entry}' not found in {fasta_path}",
              file=sys.stderr)
        return 1

    circular = mssann.is_circular(entry, common)
    ids = SeqIdMaker(args.entry, FEATURE)
    counts = {"out": 0, "total": 0, "warn": 0, "info": 0, "quoted": 0}

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        for fs in iter_feature_seqs(entry, seq, FEATURE, circular):
            feat = fs.feature
            counts["total"] += 1

            def warn(msg, _feat=feat):
                counts["warn"] += 1
                print(f"Warning: line {_feat.line_no}: {msg}", file=sys.stderr)

            def info(msg, _feat=feat):
                counts["info"] += 1
                print(f"Info: line {_feat.line_no}: {msg}", file=sys.stderr)

            if fs.error is not None:
                warn(f"{fs.error}; skipped")
                continue
            if feat.location_quoted:
                counts["quoted"] += 1
            for problem in feat.pseudo_problems():
                warn(problem)
            note = partial_note(fs.location)
            if note:
                info(note)
            check_rrna(feat, warn)

            seq_id = ids.make(feat, fs.index)
            out_seq = fs.seq
            if args.rna:
                out_seq = out_seq.replace("T", "U").replace("t", "u")

            header = (f">{seq_id} entry={args.entry} feature={feat.key} "
                      f"location={feat.location} length={len(out_seq)}")
            for key in ("gene", "product"):
                if feat.get(key):
                    header += f" {key}={feat.get(key)}"
            if feat.get("pseudogene") is not None:
                header += f" pseudogene={feat.get('pseudogene')}"
            elif feat.has("pseudo"):
                header += " pseudo"
            out.write(header + "\n" + wrap(out_seq, args.width) + "\n")
            counts["out"] += 1
    finally:
        if out is not sys.stdout:
            out.close()

    # File-level notes (reported once)
    if ids.n_without_locus_tag:
        print(f"Info: {ids.n_without_locus_tag} {FEATURE} feature(s) have no "
              "locus_tag; gene names are used as sequence IDs", file=sys.stderr)
        counts["info"] += 1
    if counts["quoted"]:
        print(f"Info: {counts['quoted']} location(s) are written in double "
              "quotes (spreadsheet quoting); quotes are removed",
              file=sys.stderr)
        counts["info"] += 1

    print(f"# {counts['out']}/{counts['total']} {FEATURE} sequence(s) from "
          f"entry '{args.entry}', {counts['warn']} warning(s), "
          f"{counts['info']} info", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
