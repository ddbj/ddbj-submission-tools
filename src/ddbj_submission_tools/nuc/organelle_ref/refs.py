"""Reference genomes: names, the NCBI search (genus only unless asked), one genome per species, cache.

Only the genus name (and, with --send-species, the species name) goes to NCBI; sequences never do.
"""
from __future__ import annotations

import http.client
import json
import os
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "ddbj-submission-tools"
_UNRESOLVED_FIRST = {"uncultured", "unidentified", "unclassified"}
_PAGE = 500


class ReferenceError(RuntimeError):
    pass


@dataclass(frozen=True)
class RefSummary:
    accession: str
    species: str
    taxid: int
    title: str

    @property
    def refseq(self) -> bool:
        return self.accession.startswith("NC_")


def genus_of(organism: str) -> str | None:
    words = organism.split()
    if not words or words[0].lower() in _UNRESOLVED_FIRST:
        return None
    return words[0]


def species_of(organism: str) -> str | None:
    """The binomial, or None when the species is not resolved (sp., cf., aff., hybrids, uncultured ...)."""
    words = organism.split()
    if len(words) < 2 or genus_of(organism) is None:
        return None
    if words[1] in ("sp.", "sp") or "cf." in words or "aff." in words or "x" in words or "×" in words:
        return None
    return f"{words[0]} {words[1]}"


def organelle_filter(value: str) -> str:
    v = value.strip().lower()
    if v.startswith("mitochondrion"):
        return "mitochondrion"
    if v == "plastid:chloroplast":
        return "chloroplast"
    if v.startswith("plastid"):
        return "plastid"
    raise ValueError(f"no complete-genome search for organelle '{value}'")


def pick_one_per_species(summaries: list[RefSummary], claimed: str | None,
                         max_species: int) -> tuple[list[RefSummary], int]:
    """One genome per species (RefSeq first), species with RefSeq first, then search order; claimed always kept."""
    best: dict[str, RefSummary] = {}
    for s in summaries:
        cur = best.get(s.species)
        if cur is None or (s.refseq and not cur.refseq):
            best[s.species] = s
    ordered = sorted(best.values(), key=lambda s: not s.refseq)          # stable: keeps search order otherwise
    picked = ordered[:max_species]
    if claimed in best and best[claimed] not in picked:
        picked = picked[:max_species - 1] + [best[claimed]] if max_species > 0 else [best[claimed]]
    return picked, len(ordered) - len(picked)


def _default_get(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read()


class NcbiClient:
    def __init__(self, get: Callable[[str], bytes] | None = None, pause: float = 0.34):
        self._get = get or _default_get
        self.pause = pause

    def _call(self, util: str, **params) -> bytes:
        params["tool"] = TOOL
        url = f"{EUTILS}/{util}.fcgi?{urllib.parse.urlencode(params)}"
        try:
            data = self._get(url)
        except (urllib.error.URLError, http.client.HTTPException, OSError) as e:
            raise ReferenceError(f"NCBI E-utilities not reachable: {e}") from e
        finally:
            time.sleep(self.pause)
        return data

    def _json(self, util: str, **params) -> dict:
        try:
            return json.loads(self._call(util, retmode="json", **params))
        except ValueError as e:
            raise ReferenceError(f"NCBI {util}: unexpected answer") from e

    def _search_term(self, term: str) -> list[RefSummary]:
        try:
            es = self._json("esearch", db="nuccore", term=term, usehistory="y", retmax=0)["esearchresult"]
            count, webenv, key = int(es["count"]), es.get("webenv", ""), es.get("querykey", "")
            out = []
            for start in range(0, count, _PAGE):
                res = self._json("esummary", db="nuccore", WebEnv=webenv, query_key=key,
                                 retstart=start, retmax=_PAGE)["result"]
                for uid in res.get("uids", []):
                    d = res[uid]
                    species = " ".join(d.get("organism", "").split()[:2])
                    out.append(RefSummary(d["accessionversion"].split(".")[0], species,
                                          int(d.get("taxid") or 0), d.get("title", "")))
            return out
        except (KeyError, TypeError, ValueError) as e:
            raise ReferenceError(f"NCBI search: unexpected answer ({e})") from e

    def search(self, genus: str, organelle: str, species: str | None = None) -> list[RefSummary]:
        tail = f" AND {organelle}[filter] AND complete genome[title]"
        found = self._search_term(f"{genus}[Organism]{tail}")
        if species:
            seen = {s.accession for s in found}
            found += [s for s in self._search_term(f"{species}[Organism]{tail}") if s.accession not in seen]
        return found

    def fetch(self, accessions: list[str], cache_dir: Path) -> list[Path]:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for acc in accessions:
            p = cache_dir / f"{acc}.gb"
            if not p.exists():
                text = self._call("efetch", db="nuccore", id=acc, rettype="gbwithparts", retmode="text")
                fd, tmp = tempfile.mkstemp(dir=cache_dir, suffix=".tmp")
                with os.fdopen(fd, "wb") as fh:
                    fh.write(text)
                os.replace(tmp, p)
            paths.append(p)
        return paths


def load_refs(paths: list[Path]) -> dict:
    """accession -> (RefSummary, SeqRecord) for GenBank files (one or more records each)."""
    from Bio import SeqIO
    out = {}
    for path in paths:
        for rec in SeqIO.parse(str(path), "genbank"):
            acc = rec.name if rec.name and rec.name != "<unknown name>" else rec.id.split(".")[0]
            taxid = 0
            for f in rec.features:
                if f.type == "source":
                    for x in f.qualifiers.get("db_xref", []):
                        if x.startswith("taxon:") and x[6:].isdigit():
                            taxid = int(x[6:])
                    break
            species = " ".join(rec.annotations.get("organism", "").split()[:2])
            out[acc] = (RefSummary(acc, species, taxid, rec.description), rec)
    return out
