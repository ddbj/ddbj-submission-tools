"""Sequences shorter than a minimum length: what removing them would lose, and removing them.

The length depends on the submission type: 1,000 bp for genome sequences (ddbj-validator ANN5045,
checked for WGS of eukaryotes and prokaryotes), 200 bp for other submissions. Shorter ones are accepted when there is a reason.

Whether to remove them is the submitter's decision; the report gives what they need to decide.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

from .annfile import AnnFile

GENOME_MIN_LENGTH = 1000


def _records(fasta: str | Path) -> list[tuple[str, list[str], int]]:
    """(name, raw lines incl. the header and a closing '//', sequence length) per record."""
    recs: list[tuple[str, list[str], int]] = []
    with open(fasta, newline="", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                recs.append((line[1:].split()[0] if line[1:].split() else "", [line], 0))
            elif recs:
                name, lines, n = recs[-1]
                lines.append(line)
                if not line.startswith("//"):
                    n += len(line.strip())
                recs[-1] = (name, lines, n)
    return recs


def short_report(ann: AnnFile, fasta: str | Path, min_length: int = GENOME_MIN_LENGTH) -> dict:
    recs = _records(fasta)
    short = [(name, n) for name, _, n in recs if n < min_length]
    names = {name for name, _ in short}
    feats = [f for f in ann.features() if f.entry != "COMMON" and f.key != "source"]
    lost = Counter(f.key for f in feats if f.entry in names)
    trna_all = Counter(f.get(ann.rows, "product") for f in feats if f.key == "tRNA")
    trna_lost = Counter(f.get(ann.rows, "product") for f in feats if f.key == "tRNA" and f.entry in names)
    total = sum(n for _, _, n in recs)
    return {
        "min_length": min_length,
        "sequences": len(recs),
        "count": len(short),
        "short": [{"entry": name, "length": n} for name, n in short],
        "short_bases": sum(n for _, n in short),
        "total_bases": total,
        "lost_features": dict(lost.most_common()),
        # products whose every copy sits on a short sequence: removing them changes the genome's tRNA set
        "lost_trna_products": sorted(p for p in trna_lost if p and trna_lost[p] == trna_all[p]),
    }


def drop_short(ann: AnnFile, fasta: str | Path, min_length: int = GENOME_MIN_LENGTH) -> str:
    """Remove short sequences: their entries from `ann` (in place) and their records from the FASTA,
    which is returned as text."""
    recs = _records(fasta)
    names = {name for name, _, n in recs if n < min_length}
    keep, entry = [], ""
    for r in ann.rows:
        if r.cols[0]:
            entry = r.cols[0]
        if entry not in names:
            keep.append(r)
    ann.rows[:] = keep
    return "".join("".join(lines) for name, lines, n in recs if name not in names)
