"""Undo a table-like .ann layout where entry and feature names are repeated on every row.

The validator reports this as a duplicate entry (ANN0120, e.g. COMMON on every row) or a
duplicate single-use feature (ANN0645, e.g. KEYWORD written once per keyword). Only what
is unambiguous is changed; everything else is returned as a review item.
"""
from __future__ import annotations

from .annfile import AnnFeature, AnnFile, AnnRow

# features an entry may hold once (ddbj-validator definitions.json, singleton_feature)
SINGLE_USE = {"SUBMITTER", "ST_COMMENT", "DIVISION", "TOPOLOGY", "KEYWORD", "DATE", "DATATYPE",
              "PRIMARY_CONTIG", "SUBCONTACT", "DBLINK"}
# REFERENCE qualifiers written once per reference (singleton_qualifiers, plus status)
_REF_SINGLE = {"title", "status", "journal", "volume", "issue", "start_page", "end_page", "pb_type",
               "ui", "media", "doi", "year"}
_HEADER = ["entry", "feature", "location"]


def collapse_repeats(ann: AnnFile) -> tuple[list[str], list[str]]:
    """Edit `ann` in place. Returns (log, review): what changed, and what a person must decide
    (line numbers refer to the edited file)."""
    log: list[str] = []
    if ann.rows and [c.strip().lower() for c in ann.rows[0].cols[:3]] == _HEADER:
        del ann.rows[0]
        log.append("removed the header row (Entry / Feature / Location ...)")

    entry = key = location = ""
    seen: dict[str, AnnRow] = {}
    separated: list[tuple[str, AnnRow]] = []
    cleared_entries = cleared_keys = 0
    for r in ann.rows:
        name = r.cols[0]
        if name and name == entry:
            r.cols[0] = ""
            cleared_entries += 1
        elif name:
            if name in seen:
                separated.append((name, r))
            seen.setdefault(name, r)
            entry, key, location = name, "", ""
        if r.cols[1]:
            if r.cols[1] == key and key in SINGLE_USE and not r.cols[0] and not r.cols[2] and not location:
                r.cols[1] = ""
                cleared_keys += 1
            else:
                key, location = r.cols[1], r.cols[2]
    if cleared_entries:
        log.append(f"cleared {cleared_entries} repeated entry name(s)")
    if cleared_keys:
        log.append(f"joined {cleared_keys} repeated single-use feature row(s) to the feature above")

    pending = _references(ann, log)
    line = {id(r): i + 1 for i, r in enumerate(ann.rows)}
    review = [f"entry {name} appears again at line {line[id(r)]} after other entries: move its rows "
              f"to the first {name} by hand" for name, r in separated]
    review += [f"REFERENCE at lines {line[id(a)]} and {line[id(b)]}: {why}; decide with "
               f"--at {line[id(b)]}=merge-previous, --at {line[id(b)]}=drop, or keep both" for a, b, why in pending]
    return log, review


def _complete(q: list[tuple[str, str]]) -> bool:
    keys = {k for k, _ in q}
    return {"title", "status"} <= keys and bool(keys & {"ab_name", "consrtm"})


def _references(ann: AnnFile, log: list[str]) -> list[tuple[AnnRow, AnnRow, str]]:
    """Merge or drop consecutive REFERENCE features where the content decides it."""
    pending = []
    merged = dropped = 0
    acc: AnnFeature | None = None
    i = 0
    while True:
        feats = ann.features()
        if i >= len(feats):
            break
        f = feats[i]
        prev = feats[i - 1] if i else None
        if f.key != "REFERENCE" or acc is None or prev is None or prev.key != "REFERENCE" \
                or prev.entry != f.entry or ann.rows[f.first].cols[0]:
            acc = f if f.key == "REFERENCE" else None
            i += 1
            continue
        a, b = acc.qualifiers(ann.rows), f.qualifiers(ann.rows)
        clash = {k for k, _ in b if k in _REF_SINGLE} & {k for k, _ in a}
        if sorted(a) == sorted(b):
            ann.delete_feature(f)
            dropped += 1
            continue                                  # same index now holds the next feature
        if not clash and not (_complete(a) and _complete(b)) and not f.location:
            ann.rows[f.first].cols[1] = ""
            merged += 1
            acc = ann.features()[i - 1]
            continue
        if clash and not (_complete(a) and _complete(b)) \
                or _complete(a) and _complete(b) and _same_title(a, b):
            pending.append((ann.rows[acc.first], ann.rows[f.first], "cannot tell one split reference from two"))
        acc = f
        i += 1
    if merged:
        log.append(f"joined {merged} REFERENCE part(s) to the reference above")
    if dropped:
        log.append(f"removed {dropped} REFERENCE identical to the one above")
    return pending


def _same_title(a, b) -> bool:
    norm = lambda q: " ".join(v for k, v in q if k == "title").lower().strip(" .")
    return norm(a) == norm(b)


def merge_previous(ann: AnnFile, line: int) -> str:
    """Join the feature starting at `line` (1-based) to the feature right above it (same entry and key)."""
    feats = ann.features()
    idx = next((n for n, f in enumerate(feats) if f.first == line - 1), None)
    if idx is None:
        raise ValueError(f"line {line} is not the first row of a feature")
    f = feats[idx]
    prev = feats[idx - 1] if idx else None
    if prev is None or prev.key != f.key or prev.entry != f.entry or ann.rows[f.first].cols[0]:
        raise ValueError(f"line {line} ({f.key}) has no {f.key} of the same entry right above it")
    if f.location and f.location != prev.location:
        raise ValueError(f"line {line} ({f.key}) has its own location {f.location}")
    ann.rows[f.first].cols[1] = ann.rows[f.first].cols[2] = ""
    return f"line {line}: joined {f.key} to the {f.key} above"
