"""Drop exon, intron, 5'UTR and 3'UTR features that can be derived from mRNA and CDS.

Converting a GFF mechanically often writes every exon, intron and UTR of a eukaryotic
genome annotation. They are not wrong, but they are not required: the mRNA (exons and
introns) and the CDS inside it (UTRs) already describe them. This removes them.

Only on a nuclear genome of a cellular organism, and only what is plainly derivable.
A feature is dropped when
  - it lies within an mRNA on the same strand in the same entry (a UTR: an mRNA that
    also has a CDS within it), and
  - it carries no qualifier beyond gene, gene_synonym, locus_tag, number and inference
    (a note or other qualifier means someone wrote it on purpose).
Skipped entirely:
  - entries whose source (or the COMMON source) has organelle, an organism name with
    virus / phage / viroid, or a COMMON DIVISION of VRL / PHG: viruses and organelles
    often have unusual transcripts.
Kept and reported: features outside any mRNA (e.g. exons of a ncRNA gene), with other
qualifiers, or on mixed strands.

Do not use it when the parts were written on purpose: a submission of single genes or
gene clusters, or genes with unusual structure. Ask the submitter before applying.

Usage:
    python -m ddbj_submission_tools.nuc.fix.derivable_features FILE.ann --dry-run
    ... --outdir DIR | --in-place   (--in-place keeps FILE.ann.orig)
    --keys exon,intron                only these feature keys (default exon,intron,5'UTR,3'UTR)
    --allow-qualifier note            also drop features with this qualifier (after checking its values)
    --all                             drop every feature of the keys, derivable or not (organelle and
                                      virus entries are still skipped); only when the submitter wants
                                      none of these features

The decisions are written to stdout as TSV (entry, feature, location, action, reason).
Exit code: 0 done, 2 usage or file error.
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

from ddbj_submission_tools.nuc.ann_fix.annfile import AnnFile
from ddbj_submission_tools.nuc.mssann.location import LocationError, parse_location

DEFAULT_KEYS = {"exon", "intron", "5'UTR", "3'UTR"}
UTR_KEYS = {"5'UTR", "3'UTR"}
PLAIN_QUALIFIERS = {"gene", "gene_synonym", "locus_tag", "number", "inference"}
_VIRUS = re.compile(r"\b(virus|phage|viroid)", re.I)


@dataclass
class Result:
    entry: str
    feature: str
    location: str
    dropped: bool
    reason: str


def _span(text: str):
    """(strand, [(start, end), ...] sorted) or None for a mixed-strand or unparsable location."""
    try:
        loc = parse_location(text)
    except LocationError:
        return None
    return None if loc.strand == 0 else (loc.strand, sorted((p.start, p.end) for p in loc.spans))


def _inside(parts, regions) -> bool:
    return all(any(a <= s and e <= b for a, b in regions) for s, e in parts)


def _gaps(parts):
    return [(e1 + 1, s2 - 1) for (_, e1), (s2, _) in zip(parts, parts[1:])]


def _skip_reason(quals: list[tuple[str, str]]) -> str | None:
    for k, v in quals:
        if k == "organelle":
            return f"organelle ({v})"
        if k == "organism" and _VIRUS.search(v):
            return f"virus ({v})"
        if k == "division" and v.upper() in ("VRL", "PHG"):
            return f"division {v}"
    return None


def drop(ann: AnnFile, keys=DEFAULT_KEYS, allow=(), everything: bool = False) -> list[Result]:
    """Delete the derivable features in place; return a Result per feature of the given keys.

    everything=True (--all) drops every feature of the keys, derivable or not; organelle and
    virus entries are still skipped.
    """
    feats = ann.features()
    common = [q for f in feats if f.entry == "COMMON" and f.key in ("source", "DIVISION")
              for q in f.qualifiers(ann.rows)]
    by_entry: dict[str, list] = {}
    for f in feats:
        by_entry.setdefault(f.entry, []).append(f)

    results, doomed = [], []
    for entry, efeats in by_entry.items():
        if entry == "COMMON":
            continue
        source = [q for f in efeats if f.key == "source" for q in f.qualifiers(ann.rows)]
        if _skip_reason(common + source):
            continue
        spans = {id(f): _span(f.location) for f in efeats if f.key in ("mRNA", "CDS")}
        cds = [s for f in efeats if f.key == "CDS" and (s := spans[id(f)])]
        mrnas = []
        for f in efeats:
            if f.key == "mRNA" and (m := spans[id(f)]):
                has_cds = any(c[0] == m[0] and _inside(c[1], m[1]) for c in cds)
                mrnas.append((m, has_cds))
        for f in efeats:
            if f.key not in keys:
                continue
            s = _span(f.location)
            extra = sorted({k for k, _ in f.qualifiers(ann.rows)} - PLAIN_QUALIFIERS - set(allow))
            if everything:
                results.append(Result(entry, f.key, f.location, True, "--all"))
                doomed.append(f)
                continue
            if s is None:
                reason = "mixed-strand or unreadable location"
            elif extra:
                reason = "has " + ", ".join(extra)
            elif not (covering := [c for m, c in mrnas if m[0] == s[0] and
                                   _inside(s[1], _gaps(m[1]) if f.key == "intron" else m[1])]):
                reason = "not within an mRNA " + ("intron" if f.key == "intron" else "exon") + " on the same strand"
            elif f.key in UTR_KEYS and not any(covering):
                reason = "the mRNA covering it has no CDS"
            else:
                results.append(Result(entry, f.key, f.location, True, "derivable from mRNA/CDS"))
                doomed.append(f)
                continue
            results.append(Result(entry, f.key, f.location, False, reason))
    for f in sorted(doomed, key=lambda f: f.first, reverse=True):
        ann.delete_feature(f)
    return results


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ann", type=Path, help="annotation file (.ann)")
    ap.add_argument("--keys", default=",".join(sorted(DEFAULT_KEYS)),
                    help="feature keys to drop, comma-separated (default exon,intron,5'UTR,3'UTR)")
    ap.add_argument("--allow-qualifier", action="append", default=[], metavar="KEY",
                    help="also drop features carrying this qualifier (e.g. note, when every note is a "
                         "converter's 'source:AUGUSTUS'); repeatable")
    ap.add_argument("--all", action="store_true",
                    help="drop every feature of the keys, also those that are not derivable or carry "
                         "other qualifiers (organelle and virus entries are still skipped)")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="report only; write nothing")
    mode.add_argument("--outdir", type=Path, help="write the new .ann to this directory")
    mode.add_argument("--in-place", action="store_true", help="overwrite the .ann, keeping FILE.ann.orig")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    keys = {k.strip() for k in args.keys.split(",") if k.strip()}
    if not keys <= DEFAULT_KEYS:
        print(f"Error: --keys takes only {', '.join(sorted(DEFAULT_KEYS))}", file=sys.stderr)
        return 2
    try:
        ann = AnnFile.read(args.ann)
    except (OSError, UnicodeDecodeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2
    if args.outdir and (args.outdir / args.ann.name).resolve() == args.ann.resolve():
        print(f"Error: --outdir is the input directory: {args.outdir}", file=sys.stderr)
        return 2
    results = drop(ann, keys, args.allow_qualifier, args.all)
    print("#entry\tfeature\tlocation\taction\treason")
    for r in results:
        print(f"{r.entry}\t{r.feature}\t{r.location}\t{'drop' if r.dropped else 'keep'}\t{r.reason}")
    n = sum(r.dropped for r in results)
    verb = "would drop" if args.dry_run else "dropped"
    print(f"# {verb} {n} feature(s); kept {len(results) - n}", file=sys.stderr)
    if n and not args.dry_run:
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
