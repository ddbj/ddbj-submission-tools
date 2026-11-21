"""Cut feature sequences out of an entry sequence.


Shared by the sequence extraction scripts (nuc.check.rrna, nuc.check.trna, nuc.inspect.cds_seq): parses each feature's location,
extracts its nucleotide sequence and assigns a sequence ID.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterator, Optional

from .annotation import Entry, Feature
from .location import Location, LocationError, parse_location


@dataclass
class FeatureSeq:
    """One feature with its parsed location and extracted sequence.

    If the location is missing or invalid, location and seq are None and
    error holds the reason.
    """
    feature: Feature
    index: int                       # 1-based index among features of the key
    location: Optional[Location]
    seq: Optional[str]
    error: Optional[str] = None


def iter_feature_seqs(entry: Entry, seq: str, key: str,
                      circular: bool = False) -> Iterator[FeatureSeq]:
    """Yield FeatureSeq for every feature of the given key in file order."""
    for idx, feat in enumerate(entry.get_features(key), start=1):
        if feat.location is None:
            yield FeatureSeq(feat, idx, None, None, "no location")
            continue
        try:
            loc = parse_location(feat.location, seq_length=len(seq),
                                 circular=circular)
            yield FeatureSeq(feat, idx, loc, loc.extract(seq))
        except LocationError as exc:
            yield FeatureSeq(feat, idx, None, None, str(exc))


class SeqIdMaker:
    """Make sequence IDs: locus_tag, else gene (made unique), else index.

    A gene that appears more than once (e.g. in the inverted repeats of a
    plastid genome) gets suffixes _2, _3, ...
    """

    def __init__(self, entry_name: str, key: str):
        self.entry_name = entry_name
        self.key = key
        self.used: Dict[str, int] = {}
        self.n_without_locus_tag = 0

    def make(self, feat: Feature, index: int) -> str:
        locus_tag = feat.get("locus_tag")
        if locus_tag is not None:
            return locus_tag
        self.n_without_locus_tag += 1
        base = feat.get("gene") or f"{self.entry_name}_{self.key}{index}"
        self.used[base] = self.used.get(base, 0) + 1
        n = self.used[base]
        return base if n == 1 else f"{base}_{n}"


def partial_note(loc: Location) -> str:
    """Return "'<' is found", "'>' is found", "'<' and '>' are found",
    or "" for a location without partial markers."""
    has_lt = any(s.partial_start for s in loc.spans)
    has_gt = any(s.partial_end for s in loc.spans)
    if has_lt and has_gt:
        return "'<' and '>' are found"
    if has_lt:
        return "'<' is found"
    if has_gt:
        return "'>' is found"
    return ""
