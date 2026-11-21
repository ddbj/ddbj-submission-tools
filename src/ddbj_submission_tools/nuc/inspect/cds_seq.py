"""Extract nucleotide (or translated amino acid) sequences of CDS features.

For each CDS of the given entry, the location is parsed with
mssann.parse_location and the sequence is cut out of the .fasta entry
with Location.extract (join / complement / complement(join) handled).

Output is FASTA. Header format:
    >ID entry=ENTRY feature=CDS location=LOCATION length=N
        [gene=...] [product=...] [pseudogene=VALUE | pseudo]
ID is the locus_tag. Without locus_tag, the gene name is used (with a
suffix _2, _3, ... when a gene appears more than once, e.g. in the
inverted repeats of a plastid genome); without gene, ENTRY_CDSn is used.

Pseudo CDS features (pseudogene qualifier, or the old pseudo qualifier)
are output as nucleotide sequences but not translated. The old "pseudo"
qualifier cannot be used in new submissions and is reported as Warning;
pseudogene must be one of: processed, unprocessed, unitary, allelic,
unknown (other values are reported as Warning).

Locations written in spreadsheet quoting ("...") are accepted; the quotes
are removed by the parser and reported once as Info.

With --translate, the sequence is translated with Biopython using the
CDS qualifiers codon_start and transl_table (default 1). Biopython is
needed only for --translate. In this mode:
  - partial locations ('<' or '>'), a missing transl_table (table 1 used),
    a length that is not a multiple of 3, and an alternative start codon
    read as M are reported as "Info:" lines;
  - an invalid transl_table or codon_start is reported as "Warning:" and
    that CDS is skipped;
  - when the 5' end is complete (no partial marker at the 5' end) and
    codon_start is 1, a start codon of the genetic code is translated as M.

Usage:
    python -m ddbj_submission_tools.nuc.inspect.cds_seq FILE.ann ENTRY_NAME [--feature CDS]
        [--fasta FILE.fasta] [--translate] [--width 60] [-o OUT.fasta]
"""

import argparse
import sys
from pathlib import Path

from ddbj_submission_tools.nuc import mssann
from ddbj_submission_tools.nuc.mssann import fasta


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract CDS sequences of an entry as FASTA.")
    parser.add_argument("ann", type=Path, help="annotation file (.ann)")
    parser.add_argument("entry", help="entry name (e.g. chr01)")
    parser.add_argument("--feature", default="CDS",
                        help="feature key to extract (default: CDS)")
    parser.add_argument("--fasta", type=Path, default=None,
                        help="sequence file (default: same base name as .ann)")
    parser.add_argument("--translate", action="store_true",
                        help="output amino acid sequences (needs Biopython)")
    parser.add_argument("--width", type=int, default=60,
                        help="line width of output sequences (0: no wrap)")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="output file (default: stdout)")
    return parser.parse_args(argv)


def wrap(seq, width):
    if width <= 0:
        return seq
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def partial_markers(loc):
    """Return a note such as "'>' is found" for partial locations, or ""."""
    has_lt = any(s.partial_start for s in loc.spans)
    has_gt = any(s.partial_end for s in loc.spans)
    if has_lt and has_gt:
        return "'<' and '>' are found"
    if has_lt:
        return "'<' is found"
    if has_gt:
        return "'>' is found"
    return ""


def translate(nt_seq, feat, loc, warn, info):
    """Translate a CDS with Biopython, honouring codon_start/transl_table.

    Returns the amino acid sequence, or None if the CDS cannot be
    translated (invalid transl_table or codon_start; reported as Warning).

    Info messages (one line, joined by "; "):
      - partial markers '<' / '>'
      - transl_table not given (table 1 is used)
      - length not a multiple of 3
      - first codon translated as M (alternative start codon)
    """
    # Imported here: only --translate needs Biopython
    from Bio.Data import CodonTable
    from Bio.Seq import Seq

    notes = []
    marker = partial_markers(loc)
    if marker:
        notes.append(marker)

    # Genetic code: the CDS's own transl_table, default 1 (Standard)
    table_value = feat.get("transl_table")
    if table_value is None:
        table = 1
        notes.append("transl_table is not found; table 1 is used")
    else:
        try:
            table = int(table_value)
        except ValueError:
            table = None
        if table not in CodonTable.unambiguous_dna_by_id:
            warn(f"invalid transl_table '{table_value}'; not translated")
            return None

    codon_start_value = feat.get("codon_start", "1")
    if codon_start_value not in ("1", "2", "3"):
        warn(f"invalid codon_start '{codon_start_value}'; not translated")
        return None
    codon_start = int(codon_start_value)

    if feat.has("transl_except"):
        warn(f"transl_except is not applied ({feat.get('transl_except')})")

    coding = nt_seq[codon_start - 1:]
    remainder = len(coding) % 3
    if remainder:
        notes.append(f"length is not a multiple of 3 ({len(coding)} bp); "
                     f"trailing {remainder} base(s) ignored")
        coding = coding[:len(coding) - remainder]

    protein = str(Seq(coding).translate(table=table))

    # A complete 5' end starts with an initiation codon, which is always
    # read as Met (same as Biopython's cds=True and DDBJ translation).
    if protein and not loc.five_prime_partial and codon_start == 1:
        first_codon = coding[:3].upper().replace("U", "T")
        start_codons = CodonTable.unambiguous_dna_by_id[table].start_codons
        if first_codon in start_codons and protein[0] != "M":
            notes.append(f"start codon {first_codon} is translated as M "
                         f"(table {table})")
            protein = "M" + protein[1:]

    if notes:
        info("; ".join(notes))

    if "*" in protein[:-1]:
        warn(f"internal stop codon(s) at "
             f"{[i + 1 for i, aa in enumerate(protein[:-1]) if aa == '*']}")
    return protein


def main(argv=None):
    args = parse_args(argv)

    if args.translate:
        try:
            import Bio  # noqa: F401  (check availability before work)
        except ImportError:
            print("Error: --translate needs Biopython "
                  "(pip install biopython)", file=sys.stderr)
            return 1

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

    # Read only the requested entry from the sequence file
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
    features = entry.get_features(args.feature)

    out = open(args.output, "w") if args.output else sys.stdout
    n_out = n_warn = n_info = 0
    n_quoted = n_no_tag = 0
    used_ids = {}
    try:
        for idx, feat in enumerate(features, start=1):
            def warn(msg, _feat=feat):
                nonlocal n_warn
                n_warn += 1
                print(f"Warning: line {_feat.line_no}: {msg}", file=sys.stderr)

            def info(msg, _feat=feat):
                nonlocal n_info
                n_info += 1
                print(f"Info: line {_feat.line_no}: {msg}", file=sys.stderr)

            if feat.location is None:
                warn("no location; skipped")
                continue
            try:
                loc = mssann.parse_location(feat.location,
                                            seq_length=len(seq),
                                            circular=circular)
                nt_seq = loc.extract(seq)
            except mssann.LocationError as exc:
                warn(f"{exc}; skipped")
                continue

            if feat.location_quoted:
                n_quoted += 1

            # Sequence ID: locus_tag, else gene (made unique), else index
            gene = feat.get("gene")
            seq_id = feat.get("locus_tag")
            if seq_id is None:
                n_no_tag += 1
                base = gene or f"{args.entry}_{args.feature}{idx}"
                used_ids[base] = used_ids.get(base, 0) + 1
                seq_id = base if used_ids[base] == 1 else f"{base}_{used_ids[base]}"

            # pseudo / pseudogene rules (checked in both output modes)
            is_pseudo = feat.is_pseudo
            for problem in feat.pseudo_problems():
                warn(problem)
            if args.translate:
                if is_pseudo:
                    info("pseudogene; not translated")
                    continue
                out_seq = translate(nt_seq, feat, loc, warn, info)
                if out_seq is None:
                    continue
            else:
                out_seq = nt_seq

            header = (f">{seq_id} entry={args.entry} feature={feat.key} "
                      f"location={feat.location} length={len(out_seq)}")
            if gene:
                header += f" gene={gene}"
            product = feat.get("product")
            if product:
                header += f" product={product}"
            pseudogene = feat.get("pseudogene")
            if pseudogene is not None:
                header += f" pseudogene={pseudogene}"
            elif feat.has("pseudo"):
                header += " pseudo"
            out.write(header + "\n" + wrap(out_seq, args.width) + "\n")
            n_out += 1
    finally:
        if out is not sys.stdout:
            out.close()

    # File-level notes (reported once, not per feature)
    if n_no_tag:
        print(f"Info: {n_no_tag} {args.feature} feature(s) have no locus_tag; "
              f"gene names are used as sequence IDs", file=sys.stderr)
        n_info += 1
    if n_quoted:
        print(f"Info: {n_quoted} location(s) are written in double quotes "
              f"(spreadsheet quoting); quotes are removed", file=sys.stderr)
        n_info += 1

    kind = "amino acid" if args.translate else "nucleotide"
    print(f"# {n_out}/{len(features)} {args.feature} {kind} sequence(s) "
          f"from entry '{args.entry}', {n_warn} warning(s), {n_info} info",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
