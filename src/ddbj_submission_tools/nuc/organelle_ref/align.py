"""Local alignment with BLAST+ (blastn) and mapping positions between a query and a reference."""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

_FMT = "6 qseqid sseqid pident bitscore qstart qend sstart send qlen qseq sseq"


def blast_available() -> bool:
    return bool(shutil.which("blastn") and shutil.which("makeblastdb"))


@dataclass
class Hsp:
    query: str
    ref: str
    identity: float
    bitscore: float
    qstart: int
    qend: int
    sstart: int
    send: int
    qlen: int
    qseq: str
    sseq: str

    @property
    def strand(self) -> int:
        return 1 if self.sstart <= self.send else -1

    @property
    def coverage(self) -> float:
        return (self.qend - self.qstart + 1) / self.qlen

    def pairs(self):
        """(query position, reference position) per alignment column; None where that side has a gap."""
        q, s = self.qstart, self.sstart
        for a, b in zip(self.qseq, self.sseq):
            yield (None if a == "-" else q), (None if b == "-" else s)
            if a != "-":
                q += 1
            if b != "-":
                s += self.strand


def q2s(h: Hsp, qpos: int) -> int | None:
    for q, s in h.pairs():
        if q == qpos:
            return s
    return None


def s2q(h: Hsp, spos: int) -> int | None:
    for q, s in h.pairs():
        if s == spos:
            return q
    return None


def project(h: Hsp, spos: int, low: bool = True) -> tuple[int, int]:
    """Query position of reference position *spos*, and how many bases it lies beyond the HSP (0 = aligned).

    Outside the HSP the reference coordinates are carried on without indels, so the position may be
    below 1 or past the query end. A reference base deleted in the query gives the first query base
    after it (*low*, the low end of a feature in query coordinates) or the last one before it.
    """
    lo, hi = sorted((h.sstart, h.send))
    if lo <= spos <= hi:
        last = h.qstart - 1
        for q, s in h.pairs():
            if q is not None:
                last = q
            if s == spos:
                return (q, 0) if q is not None else (last + 1 if low else last, 0)
    d = lo - spos if spos < lo else spos - hi
    at_qstart = (spos < lo) == (h.strand == 1)
    return (h.qstart - d if at_qstart else h.qend + d), d


def _write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n{seq}\n")


def best_hits(queries: dict[str, str], refs: dict[str, str], workdir: Path) -> dict[str, dict[str, Hsp]]:
    """query name -> reference name -> the HSP with the highest bit score (queries without hits are absent)."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    db, qf = workdir / "refs.fa", workdir / "query.fa"
    _write_fasta(db, refs)
    _write_fasta(qf, queries)
    try:
        subprocess.run(["makeblastdb", "-in", str(db), "-dbtype", "nucl", "-out", str(db)],
                       check=True, capture_output=True, text=True)
        out = subprocess.run(["blastn", "-task", "blastn", "-evalue", "1e-10", "-query", str(qf), "-db", str(db),
                              "-outfmt", _FMT], check=True, capture_output=True, text=True).stdout
    except (OSError, subprocess.CalledProcessError) as e:
        detail = getattr(e, "stderr", "") or str(e)
        raise RuntimeError(f"BLAST+ failed: {detail.strip()}") from e
    hits: dict[str, dict[str, Hsp]] = {}
    for line in out.splitlines():
        c = line.split("\t")
        h = Hsp(c[0], c[1], float(c[2]), float(c[3]), int(c[4]), int(c[5]), int(c[6]), int(c[7]), int(c[8]),
                c[9], c[10])
        cur = hits.setdefault(h.query, {}).get(h.ref)
        if cur is None or h.bitscore > cur.bitscore:
            hits[h.query][h.ref] = h
    return hits
