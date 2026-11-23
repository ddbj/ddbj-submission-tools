"""Mark the partial ends of CDS and mRNA features in an .ann file.

CDS: an end that is not a codon of the feature's genetic code becomes partial
  - 5' end: the first codon is not a start codon, or codon_start is 2 or 3 -> '<' (or '>' on complement)
  - 3' end: the last codon is not a stop codon, or the length after codon_start is not a
    multiple of 3. A CDS with transl_except (aa:TERM) is left alone: its stop is completed
    by polyadenylation.
  Start and stop codons come from transl_table (Biopython's NCBI tables; 1 when there is none).
  pseudo / pseudogene CDS and codons with ambiguous bases are reported, not changed.

mRNA: an mRNA has UTRs outside its CDS. When an end of an mRNA is the same position as the
  same end of a CDS on the same strand inside it, the mRNA has no UTR there and that end is
  partial. This is checked after the CDS, so it follows ends the CDS step made partial.

Ends that are already partial are kept. Only the Location column changes.

Usage:
    python -m ddbj_submission_tools.nuc.fix.partial_ends FILE.ann FILE.fasta --dry-run
    ... --outdir DIR | --in-place   (--in-place keeps FILE.ann.orig)

The changes are written to stdout as TSV (entry, feature, old, new, reason); new is empty
for a CDS that was reported but not changed.
Exit code: 0 done, 2 usage or file error.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import dataclass, replace
from pathlib import Path

from Bio.Data import CodonTable

from ddbj_submission_tools.nuc.ann_fix.annfile import AnnFile
from ddbj_submission_tools.nuc.common.fasta import parse_fasta_sequences
from ddbj_submission_tools.nuc.mssann.location import Location, LocationError, parse_location


@dataclass
class Change:
    entry: str
    feature: str
    old: str
    new: str | None      # None: reported, not changed
    reason: str


def render(loc: Location) -> str:
    if loc.operator is None:
        return str(loc.spans[0])
    if loc.complement:
        inner = ",".join(str(replace(s, complement=False)) for s in loc.spans)
        return f"complement({loc.operator}({inner}))"
    return f"{loc.operator}({','.join(str(s) for s in loc.spans)})"


def five_prime(loc: Location) -> int:
    s = loc.product_spans()[0]
    return s.end if s.complement else s.start


def three_prime(loc: Location) -> int:
    s = loc.product_spans()[-1]
    return s.start if s.complement else s.end


def mark(loc: Location, five: bool = False, three: bool = False) -> None:
    if five:
        s = loc.product_spans()[0]
        if s.complement:
            s.partial_end = True
        else:
            s.partial_start = True
    if three:
        s = loc.product_spans()[-1]
        if s.complement:
            s.partial_start = True
        else:
            s.partial_end = True


def _parse(text: str) -> Location | None:
    """The location, or None when it is mixed-strand, a site, or cannot be written back as is."""
    try:
        loc = parse_location(text)
    except LocationError:
        return None
    if loc.strand == 0 or loc.is_between or render(loc) != text:
        return None
    return loc


def _cds_ends(feat, rows, seq: str) -> tuple[bool, bool, list[str]]:
    """(5' partial needed, 3' partial needed, reasons) for one CDS; reasons with no flag are notes."""
    if feat.get(rows, "pseudo") is not None or feat.get(rows, "pseudogene") is not None:
        return False, False, []
    table_id = feat.get(rows, "transl_table") or "1"
    try:
        table = CodonTable.unambiguous_dna_by_id[int(table_id)]
    except (ValueError, KeyError):
        return False, False, [f"invalid transl_table '{table_id}'"]
    try:
        codon_start = int(feat.get(rows, "codon_start") or 1)
    except ValueError:
        codon_start = 1
    loc = parse_location(feat.location)
    try:
        nt = loc.extract(seq).upper()
    except LocationError as e:
        return False, False, [str(e)]
    five = three = False
    reasons = []
    if not loc.five_prime_partial:
        first = nt[:3]
        if codon_start != 1:
            five = True
            reasons.append(f"codon_start={codon_start}")
        elif set(first) - set("ACGT"):
            reasons.append(f"5' codon {first} is ambiguous")
        elif first not in table.start_codons:
            five = True
            reasons.append(f"5' codon {first} is not a start codon (transl_table {table_id})")
    excepted = any("aa:TERM" in v.replace(" ", "") for k, v in feat.qualifiers(rows) if k == "transl_except")
    if not loc.three_prime_partial and not excepted:
        last = nt[-3:]
        if (len(nt) - (codon_start - 1)) % 3:
            three = True
            reasons.append("length is not a multiple of 3 (no complete stop codon)")
        elif set(last) - set("ACGT"):
            reasons.append(f"3' codon {last} is ambiguous")
        elif last not in table.stop_codons:
            three = True
            reasons.append(f"3' codon {last} is not a stop codon (transl_table {table_id})")
    return five, three, reasons


def fix(ann: AnnFile, seqs: dict[str, str]) -> list[Change]:
    """Mark the partial ends in place and return what was changed or reported."""
    changes: list[Change] = []
    feats = ann.features()
    cds_by_entry: dict[str, list[Location]] = {}

    def set_location(f, loc, reason):
        new = render(loc)
        ann.rows[f.first].cols[2] = new
        changes.append(Change(f.entry, f.key, f.location, new, reason))
        f.location = new

    for f in feats:
        if f.key != "CDS" or (loc := _parse(f.location)) is None:
            continue
        seq = seqs.get(f.entry)
        if seq is None:
            changes.append(Change(f.entry, f.key, f.location, None, "no sequence in the FASTA"))
        else:
            five, three, reasons = _cds_ends(f, ann.rows, seq)
            if five or three:
                mark(loc, five, three)
                set_location(f, loc, "; ".join(reasons))
            elif reasons:
                changes.append(Change(f.entry, f.key, f.location, None, "; ".join(reasons)))
        cds_by_entry.setdefault(f.entry, []).append(loc)

    for f in feats:
        if f.key != "mRNA" or (m := _parse(f.location)) is None:
            continue
        five = three = False
        for c in cds_by_entry.get(f.entry, []):
            if c.strand != m.strand or c.start < m.start or c.end > m.end:
                continue
            five |= not m.five_prime_partial and five_prime(c) == five_prime(m)
            three |= not m.three_prime_partial and three_prime(c) == three_prime(m)
        if five or three:
            mark(m, five, three)
            which = " and ".join(w for w, on in (("5'", five), ("3'", three)) if on)
            set_location(f, m, f"{which} end is the end of a CDS (no UTR)")
    return changes


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ann", type=Path, help="annotation file (.ann)")
    ap.add_argument("fasta", type=Path, help="sequence file (FASTA)")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="report only; write nothing")
    mode.add_argument("--outdir", type=Path, help="write the new .ann to this directory")
    mode.add_argument("--in-place", action="store_true", help="overwrite the .ann, keeping FILE.ann.orig")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        ann = AnnFile.read(args.ann)
        seqs = parse_fasta_sequences(str(args.fasta))
    except (OSError, UnicodeDecodeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2
    if args.outdir and (args.outdir / args.ann.name).resolve() == args.ann.resolve():
        print(f"Error: --outdir is the input directory: {args.outdir}", file=sys.stderr)
        return 2
    changes = fix(ann, seqs)
    print("#entry\tfeature\told\tnew\treason")
    for c in changes:
        print(f"{c.entry}\t{c.feature}\t{c.old}\t{c.new or ''}\t{c.reason}")
    changed = sum(c.new is not None for c in changes)
    verb = "would change" if args.dry_run else "changed"
    print(f"# {verb} {changed} location(s); {len(changes) - changed} CDS reported without change",
          file=sys.stderr)
    if changed and not args.dry_run:
        if args.in_place:
            shutil.copy2(args.ann, str(args.ann) + ".orig")
            target = args.ann
        else:
            args.outdir.mkdir(parents=True, exist_ok=True)
            target = args.outdir / args.ann.name
        ann.write(target)
    return 0


if __name__ == "__main__":
    sys.exit(main())
