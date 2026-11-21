"""Extract tRNA sequences of an entry as FASTA and check their naming.

For each tRNA feature, the location is parsed and the sequence is cut out
of the .fasta entry (join with introns and complement handled; the output
is the mature, spliced sequence on the tRNA's own strand).

Output header:
    >ID entry=ENTRY feature=tRNA location=LOCATION length=N
        [gene=...] [product=...] [anticodon=NNN] [pseudogene=VALUE | pseudo]
ID is the locus_tag, else the gene name (made unique with _2, _3, ...).

Checks (stderr):
  Warning: product not in "tRNA-Xxx" form; gene amino acid and product
           disagree; anticodon qualifier malformed or disagreeing with gene
           name / product; anticodon from the gene name not found in the
           tRNA sequence; location errors; pseudo / pseudogene problems.
  Info:    gene name not in trnX-NNN form (anticodon not checked);
           missing locus_tag and quoted locations (once per file).

Usage:
    python -m ddbj_submission_tools.nuc.check.trna FILE.ann ENTRY_NAME [--fasta FILE.fasta]
        [--rna] [--width 60] [-o OUT.fasta]
"""

import argparse
import sys
from pathlib import Path

from ddbj_submission_tools.nuc import mssann
from ddbj_submission_tools.nuc.mssann import fasta, trna
from ddbj_submission_tools.nuc.mssann.seqfeature import SeqIdMaker, iter_feature_seqs

FEATURE = "tRNA"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract tRNA sequences of an entry as FASTA.")
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


def check_trna(feat, nt_seq, seq, circular, warn, info):
    """Check gene / product / anticodon of one tRNA.

    Returns the anticodon (DNA, upper case) to show in the header, or None.
    """
    gene = feat.get("gene")
    product = feat.get("product")

    gene_aa = gene_anticodon = None
    if gene is not None:
        parsed = trna.parse_trna_gene(gene)
        if parsed is None:
            info(f"gene '{gene}' is not in trnX-NNN form; "
                 "anticodon not checked against the gene name")
        else:
            gene_aa, gene_anticodon = parsed

    product_aa = None
    if product is None:
        warn("no product")
    else:
        product_aa = trna.parse_trna_product(product)
        if product_aa is None:
            warn(f"product '{product}' is not in 'tRNA-Xxx' form")
        elif gene_aa is not None and gene_aa != product_aa:
            warn(f"gene '{gene}' ({gene_aa}) and product '{product}' disagree")

    # The anticodon qualifier gives the exact position on the entry
    anticodon = None
    value = feat.get("anticodon")
    if value is not None:
        parsed = trna.parse_anticodon(value)
        if parsed is None:
            warn(f"malformed anticodon '{value}'")
        else:
            pos_text, ac_aa = parsed
            try:
                ac_loc = mssann.parse_location(pos_text, seq_length=len(seq),
                                               circular=circular)
                anticodon = ac_loc.extract(seq).upper()
            except mssann.LocationError as exc:
                warn(f"anticodon position: {exc}")
            if anticodon is not None and len(anticodon) != 3:
                warn(f"anticodon position covers {len(anticodon)} bases")
            if product_aa is not None and ac_aa != product_aa:
                warn(f"anticodon aa:{ac_aa} and product '{product}' disagree")
            if (anticodon and gene_anticodon
                    and anticodon != gene_anticodon):
                warn(f"anticodon at {pos_text} is {anticodon} but gene "
                     f"'{gene}' says {gene_anticodon}")
    elif gene_anticodon is not None:
        # No position given: at least the anticodon must occur in the tRNA
        if gene_anticodon not in nt_seq.upper():
            warn(f"anticodon {gene_anticodon} (from gene '{gene}') is not "
                 "found in the tRNA sequence")
        anticodon = gene_anticodon

    return anticodon


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

            seq_id = ids.make(feat, fs.index)
            anticodon = check_trna(feat, fs.seq, seq, circular, warn, info)

            out_seq = fs.seq
            if args.rna:
                out_seq = out_seq.replace("T", "U").replace("t", "u")

            header = (f">{seq_id} entry={args.entry} feature={feat.key} "
                      f"location={feat.location} length={len(out_seq)}")
            for key in ("gene", "product"):
                if feat.get(key):
                    header += f" {key}={feat.get(key)}"
            if anticodon:
                if args.rna:
                    anticodon = anticodon.replace("T", "U")
                header += f" anticodon={anticodon}"
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
