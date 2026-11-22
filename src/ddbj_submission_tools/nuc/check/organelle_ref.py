"""Check organelle gene sequences against complete organelle genomes of the same genus.

For every entry whose source has /organelle (mitochondrion..., plastid...):

  1. Complete genomes of the genus are searched at NCBI (only the genus name is sent;
     --send-species also sends the species name). One genome per species, RefSeq first,
     at most --max-species species; the claimed species is always kept. Genomes are cached.
     --refs FILE.gb ... uses local GenBank files instead (nothing is sent).
  2. The entry is aligned to every genome locally with blastn (sequences never leave
     this machine).
  3. Identification: the identity to each species. SPECIES_MISMATCH / OTHER_SPECIES_MATCH /
     LOW_IDENTITY / PARTIAL_MAP point to a possible misidentification or contamination;
     the numbers are guides, a person decides.
  4. Boundaries: each CDS / rRNA / tRNA is compared with the same feature of the claimed
     species' genome (or the closest one): BOUNDARY, CODON_START, INCOMPLETE_STOP (info).
     From RefSeq genomes, fixes are written with --rules as nuc.fix.edit rules, all
     disabled: enable the ones you accept, then
       python -m ddbj_submission_tools.nuc.fix.edit RULES.json FILE.ann --outdir fixed

Needs BLAST+ (blastn, makeblastdb) on PATH, or run it in the container image.

Usage:
    python -m ddbj_submission_tools.nuc.check.organelle_ref FILE.ann FILE.fasta [--rules proposals.json] [--tsv]

Exit code: 0 no finding (information only), 1 findings, 2 usage, file, network or BLAST+ error.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

from ddbj_submission_tools.nuc.mssann.annotation import AnnParseError, effective_source, read
from ddbj_submission_tools.nuc.mssann.fasta import iter_sequences
from ddbj_submission_tools.nuc.organelle_ref import align, refs
from ddbj_submission_tools.nuc.organelle_ref.judge import FINDING_RESULTS, Finding, Thresholds, boundary, hit_table, identify

COLUMNS = ["entry", "check", "result", "detail", "ref", "proposal"]
DEFAULT_CACHE = Path.home() / ".cache" / "ddbj-submission-tools" / "organelle"


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ann", type=Path)
    ap.add_argument("fasta", type=Path)
    ap.add_argument("--refs", type=Path, nargs="+", help="use these GenBank files as references (no NCBI search)")
    ap.add_argument("--genus", help="search this genus instead of the first word of the organism")
    ap.add_argument("--send-species", action="store_true",
                    help="also send the species name to NCBI (default: the genus only)")
    ap.add_argument("--max-species", type=int, default=50)
    ap.add_argument("--cache", type=Path, default=DEFAULT_CACHE, help=f"GenBank cache (default {DEFAULT_CACHE})")
    ap.add_argument("--min-coverage", type=float, default=0.95)
    ap.add_argument("--delta", type=float, default=1.0)
    ap.add_argument("--other-identity", type=float, default=99.0)
    ap.add_argument("--min-identity", type=float, default=97.0)
    ap.add_argument("--tsv", action="store_true", help="write the findings as TSV to stdout (report to stderr)")
    ap.add_argument("--rules", type=Path, help="write the proposed fixes as nuc.fix.edit rules (all disabled)")
    return ap.parse_args(argv)


def _targets(ann: Path, fasta: Path):
    """[(name, seq, organism, organelle, features)] for entries with /organelle."""
    af = read(ann)
    out = []
    for name, seq in iter_sequences(fasta):
        entry = af.get_entry(name)
        src = effective_source(entry, af.common)
        organelle = src.get("organelle") if src else None
        if organelle:
            out.append((name, seq.upper(), src.get("organism") or "", organelle, entry.features if entry else []))
    return out


def _references(args, groups, client):
    """(genus, filter) -> {accession: (RefSummary, SeqRecord)}, plus information findings."""
    info: list[Finding] = []
    if args.refs:
        loaded = refs.load_refs(args.refs)
        return {g: loaded for g in groups}, info
    out = {}
    for (genus, filt), claimed in groups.items():
        if genus is None:
            out[(genus, filt)] = {}
            continue
        species = sorted(claimed)[0] if args.send_species and claimed else None
        found = client.search(genus, filt, species)
        picked, dropped = refs.pick_one_per_species(found, claimed, args.max_species)
        if dropped:
            info.append(Finding("*", "identify", "CAP_REACHED",
                                f"{genus} {filt}: {dropped} species left out by --max-species {args.max_species}; "
                                "raise it if the closest species may be among them"))
        paths = client.fetch([p.accession for p in picked], args.cache)
        out[(genus, filt)] = refs.load_refs(paths)
    return out, info


def run(args) -> tuple[list[Finding], list[dict]]:
    targets = _targets(args.ann, args.fasta)
    if not targets:
        return [], []
    if not align.blast_available():
        raise RuntimeError("blastn and makeblastdb (BLAST+) are not on PATH; install BLAST+ or use the container image")
    keyed = []
    groups: dict = {}
    for name, seq, organism, organelle, features in targets:
        filt = refs.organelle_filter(organelle)
        genus = args.genus or refs.genus_of(organism)
        claimed = refs.species_of(organism)
        groups.setdefault((genus, filt), set()).update({claimed} if claimed else set())
        keyed.append(((genus, filt), name, seq, claimed, features))
    by_group, findings = _references(args, groups, refs.NcbiClient())
    th = Thresholds(args.min_coverage, args.delta, args.other_identity, args.min_identity)
    rules: list[dict] = []
    with tempfile.TemporaryDirectory() as tmp:
        for i, (key, name, seq, claimed, features) in enumerate(keyed):
            references = by_group[key]
            if not references:
                findings.append(Finding(name, "identify", "NO_REFERENCE",
                                        f"no complete {key[1]} genome of the genus {key[0] or '(unknown)'}"))
                continue
            hits = align.best_hits({name: seq}, {acc: str(rec.seq) for acc, (_, rec) in references.items()},
                                   Path(tmp) / str(i)).get(name, {})
            ranked = sorted(((references[acc][0], h) for acc, h in hits.items()),
                            key=lambda rh: (-rh[1].identity, -rh[1].bitscore))
            findings += identify(name, claimed, ranked, th)
            if not ranked:
                continue
            full = [rh for rh in ranked if rh[1].coverage >= th.min_coverage] or ranked
            mine = [rh for rh in full if rh[0].species == claimed]
            summary, h = (mine or full)[0]
            found, new_rules = boundary(name, features, len(seq), h, references[summary.accession])
            findings += found
            rules += new_rules
            findings.append(Finding(name, "identify", "TOP", hit_table(ranked)))
    return findings, rules


def report(findings: list[Finding], out) -> None:
    current = None
    for f in findings:
        if f.entry != current:
            current = f.entry
            print(f"== {f.entry}", file=out)
        if f.result == "TOP":
            print(f"   closest: {f.detail}", file=out)
            continue
        print(f"   {f.result}\t{f.detail}" + (f"\t-> {f.proposal}" if f.proposal else ""), file=out)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        findings, rules = run(args)
    except (OSError, AnnParseError, ValueError, refs.ReferenceError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2
    out = sys.stderr if args.tsv else sys.stdout
    if not findings:
        print("no organelle entries", file=out)
    report(findings, out)
    if args.tsv:
        print("\t".join(COLUMNS))
        for f in findings:
            if f.result != "TOP":
                print("\t".join([f.entry, f.check, f.result, f.detail, f.ref, f.proposal]))
    if args.rules:
        args.rules.write_text(json.dumps({"rules": rules}, ensure_ascii=False, indent=2) + "\n")
        print(f"# {len(rules)} proposed rule(s), all disabled, in {args.rules}", file=sys.stderr)
    return 1 if any(f.result in FINDING_RESULTS for f in findings) else 0


if __name__ == "__main__":
    sys.exit(main())
