"""Write metadata back where the caller says (`ann_meta apply --placed`).

The input (``placed.json``) says where source and DBLINK go and gives the rows
as they are to be written::

    {"SUBMITTER": {...}, "REFERENCE": [...], "COMMENT": [...], "ST_COMMENT": [...],
     "_extra": {feature: [[name, value], ...]},          # optional, as in the common JSON
     "DATE": {...},                                       # optional, ignored
     "source_layout": "common" | "entry",
     "dblink_layout": "common" | "entry" | "none",
     "common_source": [[name, value], ...],               # source_layout = common
     "common_dblink": {"project": ..., "sample": ..., "sequence_read_archive": ...},  # dblink_layout = common
     "entries": [{"seq_id": ..., "topology": "circular" | "",
                  "source": [[name, value], ...],         # source_layout = entry
                  "dblink": {...}}]}                      # dblink_layout = entry

A qualifier value ``true`` or ``""`` is written as a name with an empty value
column. Rows are written exactly as given: no category, no automatic
qualifiers, no role qualifier, no ff_definition regeneration.

What is rewritten (every other row keeps its bytes):

* COMMON DBLINK (``common_dblink`` in layout ``common``, else removed),
  SUBMITTER, REFERENCE, COMMENT and ST_COMMENT, built with
  ``create_common`` (ST_COMMENT with its tags in the input's order, not
  sorted; a repeated tag from ``_extra`` right after the last row of that tag). Old and new features are matched per key by occurrence
  (whatever order the file has): a matched one is rewritten in place, an
  unmatched old one is removed (an absent or empty section), and a new one
  goes next to its neighbours in that canonical order (else before the
  COMMON source, else at the end of COMMON). DATE is never touched.
* COMMON source at ``1..E`` (layout ``common``) or removed (``entry``).
* Per entry: DBLINK (layout ``entry``, else removed), TOPOLOGY (``circular``,
  or a circular one removed; a non-circular one such as ``linear`` stays) and source at ``1..<length>`` (layout ``entry``
  with source rows; else removed: an entry given no rows in layout ``entry`` is written without a source). A new TOPOLOGY goes before the entry's source (or first in the
  entry), a new DBLINK directly before the source (or after a leading
  TOPOLOGY), a new source after the entry's leading TOPOLOGY / DBLINK. An
  entry without rows in the .ann (a FASTA entry) that needs rows gets them
  before the first entry that follows it in FASTA order and is in the .ann.

A feature whose rows already say what the input says is left as it is, so
saving without a change gives the same bytes. For source the rows and their
order must match; for ST_COMMENT the rows in order; for the other COMMON
features and DBLINK the values per qualifier name (in order), because
``create_common`` writes the names in its own order. Line endings come from the original file.
"""
from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field

from ddbj_submission_tools.nuc.ann_fix.annfile import AnnFeature, AnnFile, AnnRow
from ddbj_submission_tools.nuc.ann_meta.apply import common_blocks
from ddbj_submission_tools.nuc.ann_meta.shape import normalize_common, shape_issues
from ddbj_submission_tools.nuc.ann_meta.validate import Issue, char_issues

SOURCE_LAYOUTS = ("common", "entry")
DBLINK_LAYOUTS = ("common", "entry", "none")
TOPOLOGIES = ("circular", "")
COMMON_SECTIONS = ("SUBMITTER", "REFERENCE", "COMMENT", "ST_COMMENT")
TARGETS = ("DBLINK", *COMMON_SECTIONS)          # the order of apply's COMMON targets, without DATE
_TOP_KEYS = frozenset({*COMMON_SECTIONS, "DATE", "_extra", "source_layout", "dblink_layout",
                       "common_source", "common_dblink", "entries"})
_ENTRY_KEYS = frozenset({"seq_id", "topology", "source", "dblink"})
_QUALIFIER_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]*")

Quals = list[tuple[str, str]]


@dataclass
class PlacedEntry:
    seq_id: str
    topology: str = ""
    source: Quals = field(default_factory=list)
    dblink: dict = field(default_factory=dict)


@dataclass
class PlacedDoc:
    common: dict                    # normalized SUBMITTER / REFERENCE / COMMENT / ST_COMMENT / DBLINK / _extra
    source_layout: str
    dblink_layout: str
    common_source: Quals
    entries: list[PlacedEntry]

    @classmethod
    def from_json(cls, data: dict) -> "PlacedDoc":
        """Build from a placed JSON that :func:`placed_issues` accepted."""
        source_layout, dblink_layout = data["source_layout"], data["dblink_layout"]
        common = {k: data[k] for k in COMMON_SECTIONS if data.get(k) is not None}
        if dblink_layout == "common" and data.get("common_dblink"):
            common["DBLINK"] = data["common_dblink"]
        extra = {k: v for k, v in (data.get("_extra") or {}).items()
                 if k in COMMON_SECTIONS or (k == "DBLINK" and dblink_layout == "common")}
        if extra:
            common["_extra"] = extra
        common = normalize_common(common)
        entries = []
        for e in data["entries"]:
            entries.append(PlacedEntry(
                seq_id=e["seq_id"], topology=e.get("topology") or "",
                source=_rows(e.get("source")) if source_layout == "entry" else [],
                dblink=normalize_common({"DBLINK": e.get("dblink") or {}})["DBLINK"]
                if dblink_layout == "entry" else {}))
        return cls(common=common, source_layout=source_layout, dblink_layout=dblink_layout,
                   common_source=_rows(data.get("common_source")) if source_layout == "common" else [],
                   entries=entries)


def _rows(pairs) -> Quals:
    """[[name, value], ...] -> qualifier rows; ``true`` is a valueless qualifier."""
    return [(name, "" if value is True else value) for name, value in pairs or []]


# -- checks ------------------------------------------------------------------------

def _segments(ann: AnnFile) -> list[str]:
    """The entry marker of every run of rows that starts with one, in file order."""
    return [r.cols[0] for r in ann.rows if r.cols[0]]


def _pair_issues(value, path: str, entry: str | None, *, required: bool) -> list[Issue]:
    if value is None or value == []:
        return [Issue(path, entry, "missing_source")] if required else []
    if not isinstance(value, list):
        return [Issue(path, entry, "bad_shape: must be a list of [name, value] pairs")]
    issues = []
    for i, pair in enumerate(value):
        if not (isinstance(pair, list) and len(pair) == 2 and isinstance(pair[0], str)
                and (isinstance(pair[1], str) or pair[1] is True)):
            issues.append(Issue(f"{path}.{i}", entry,
                                "bad_shape: must be [name, value] with a string name and a string or true value"))
        elif not _QUALIFIER_NAME.fullmatch(pair[0]):
            issues.append(Issue(f"{path}.{i}", entry, "bad_qualifier_name"))
    return issues


def _dblink_issues(value, path: str, entry: str | None) -> list[Issue]:
    if value is None:
        return []
    if not isinstance(value, dict):
        return [Issue(path, entry, "bad_shape: must be an object")]
    return [Issue(path + i.field[len("DBLINK"):], entry, i.message)
            for i in shape_issues({"DBLINK": value})]


def placed_issues(data: dict, ann: AnnFile, lengths: dict[str, int],
                  fasta_ids: list[str] | None) -> list[Issue]:
    """Every problem of the placed JSON *data* for *ann* (empty: :func:`apply_placed` may run).

    *lengths* are the sequence lengths known (FASTA, else the original ``1..N``); *fasta_ids*
    the FASTA's ids when one is given. The entries must be exactly the entries of the .ann
    and of the FASTA.
    """
    issues = char_issues(data)
    if issues:
        return issues
    issues = [Issue(k, None, "not_allowed_with_placed") for k in data if k not in _TOP_KEYS]
    issues += shape_issues({k: data[k] for k in (*COMMON_SECTIONS, "_extra") if k in data})

    source_layout, dblink_layout = data.get("source_layout"), data.get("dblink_layout")
    if source_layout not in SOURCE_LAYOUTS:
        issues.append(Issue("source_layout", None, "bad_source_layout"))
    if dblink_layout not in DBLINK_LAYOUTS:
        issues.append(Issue("dblink_layout", None, "bad_dblink_layout"))
    if source_layout == "common":
        issues += _pair_issues(data.get("common_source"), "common_source", None, required=True)
    if dblink_layout == "common":
        issues += _dblink_issues(data.get("common_dblink"), "common_dblink", None)

    for m, n in Counter(_segments(ann)).items():
        if n < 2:
            continue
        issues.append(Issue("ann", None if m == "COMMON" else m, "split_entry"))

    entries = data.get("entries")
    if not isinstance(entries, list):
        return issues + [Issue("entries", None, "bad_shape: must be a list")]
    in_ann = list(dict.fromkeys(f.entry for f in ann.features() if f.entry and f.entry != "COMMON"))
    with_source = {f.entry for f in ann.features() if f.key == "source"}
    expected = dict.fromkeys(in_ann + list(fasta_ids or []))
    seen: set[str] = set()
    for i, e in enumerate(entries):
        path = f"entries.{i}"
        if not isinstance(e, dict):
            issues.append(Issue(path, None, "bad_shape: must be an object"))
            continue
        sid = e.get("seq_id")
        if not isinstance(sid, str):
            issues.append(Issue(f"{path}.seq_id", None, "bad_shape: must be a string"))
            continue
        eid = sid or None
        issues += [Issue(f"{path}.{k}", eid, "not_allowed_with_placed") for k in e if k not in _ENTRY_KEYS]
        if e.get("topology", "") not in TOPOLOGIES:
            issues.append(Issue(f"{path}.topology", eid, "bad_topology"))
        if source_layout == "entry":
            # no source rows: the entry is written without a source (COMMON's stays required)
            issues += _pair_issues(e.get("source"), f"{path}.source", eid, required=False)
            if e.get("source") and sid and sid not in with_source and sid not in lengths:
                issues.append(Issue("entries", sid, "length_unknown"))
        if dblink_layout == "entry":
            issues += _dblink_issues(e.get("dblink"), f"{path}.dblink", eid)
        if not sid:
            issues.append(Issue(f"{path}.seq_id", None, "empty_seq_id"))
        elif sid in seen:
            issues.append(Issue(f"{path}.seq_id", sid, "duplicate_seq_id"))
        elif sid not in expected:
            seen.add(sid)
            issues.append(Issue("entries", sid, "unknown_entry"))
        else:
            seen.add(sid)
    for sid in expected:
        if sid not in seen:
            issues.append(Issue("entries", sid, "missing_entry"))
    return issues


# -- writing ------------------------------------------------------------------------

class _Part(AnnFile):
    """The rows of one entry (or of COMMON); new rows take the file's line ending."""

    def __init__(self, rows: list[AnnRow], eol: str):
        super().__init__(rows)
        self._eol = eol

    def default_eol(self) -> str:
        return self._eol


def _cells(ann: AnnFile, f: AnnFeature) -> Quals:
    return [(r.cols[3], r.cols[4]) for r in ann.rows[f.first:f.last + 1]]


def _same_rows(ann: AnnFile, f: AnnFeature, location: str, quals: Quals) -> bool:
    return f.location == location and _cells(ann, f) == (quals or [("", "")])


def _by_name(quals: Quals) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for k, v in quals:
        if k or v:
            out.setdefault(k, []).append(v)
    return out


def _same_values(ann: AnnFile, f: AnnFeature, location: str, quals: Quals) -> bool:
    """Same location and the same values per qualifier name (in order). ST_COMMENT, whose tag
    order is part of what is written, must also have its rows in the same order."""
    if f.location != location:
        return False
    if f.key == "ST_COMMENT":
        return [c for c in _cells(ann, f) if any(c)] == [q for q in quals if any(q)]
    return _by_name(_cells(ann, f)) == _by_name(quals)


def _of(part: AnnFile, key: str) -> list[AnnFeature]:
    return [f for f in part.features() if f.key == key]


def _delete(part: AnnFile, feats: list[AnnFeature]) -> None:
    for f in sorted(feats, key=lambda f: f.first, reverse=True):
        part.delete_feature(f)


def _keep_one(part: AnnFile, key: str) -> AnnFeature | None:
    """Delete all but the first *key* feature; return that one."""
    feats = _of(part, key)
    _delete(part, feats[1:])
    return feats[0] if feats else None


def _common_targets(part: AnnFile) -> list[AnnFeature]:
    return [f for f in part.features() if f.key in TARGETS]


def _apply_common_targets(part: AnnFile, blocks: list[tuple[str, str, Quals]]) -> None:
    """Match the old and the new COMMON targets per key by occurrence, whatever order the file has.

    A matched feature is rewritten in place only when it differs; an old one without a match is
    deleted; a new one is inserted next to its canonical neighbours that are in the file (before
    the next one, else after the previous one), else before the COMMON source, else at the end.
    """
    old = _common_targets(part)
    by_key: dict[str, list[int]] = {}
    for i, f in enumerate(old):
        by_key.setdefault(f.key, []).append(i)
    match: dict[int, int] = {}                   # block index -> old index
    seen: dict[str, int] = {}
    for b, (key, _, _) in enumerate(blocks):
        n = seen.get(key, 0)
        if n < len(by_key.get(key, [])):
            match[b] = by_key[key][n]
        seen[key] = n + 1
    block_of = {o: b for b, o in match.items()}
    for o in reversed(range(len(old))):          # bottom up: earlier rows keep their index
        f = old[o]
        if o not in block_of:
            part.delete_feature(f)
            continue
        key, location, quals = blocks[block_of[o]]
        if not _same_values(part, f, location, quals):
            part.replace_feature(f, key, location, quals)

    present = [block_of[o] for o in sorted(block_of)]     # block index of each target, in file order
    for b, (key, location, quals) in enumerate(blocks):
        if b in match:
            continue
        feats = _common_targets(part)
        later = [p for p in present if p > b]
        earlier = [p for p in present if p < b]
        if later:
            i = present.index(min(later))
            at = feats[i].first
        elif earlier:
            i = present.index(max(earlier)) + 1
            at = feats[i - 1].last + 1
        else:
            srcs = _of(part, "source")
            i, at = len(present), (srcs[0].first if srcs else len(part.rows))
        present.insert(i, b)
        part.insert_feature(at, "COMMON" if at == 0 else "", key, location, quals)


def _st_comment_blocks(common: dict) -> list[tuple[str, str, Quals]]:
    """ST_COMMENT blocks with the tags in the input's order (``create_common`` would sort them); a
    list value is one row joined with "; " as there. ``_extra`` rows follow, as in ``common_blocks``."""
    sections = common.get("ST_COMMENT") or []
    if isinstance(sections, dict):
        sections = [sections]
    blocks: list[tuple[str, str, Quals]] = []
    for sec in sections:
        quals = [(str(k), "; ".join(str(x) for x in v) if isinstance(v, list) else str(v))
                 for k, v in sec.items()]
        if quals:
            blocks.append(("ST_COMMENT", "", quals))
    extra = [(str(k), "" if v is None or v is True else str(v))
             for k, v in (common.get("_extra") or {}).get("ST_COMMENT", []) if v is not False]
    if extra and not blocks:
        return [("ST_COMMENT", "", extra)]
    for row in extra:
        # a repeated tag goes right after the last row of that tag (where extract found it), else last
        quals = blocks[0][2]
        at = max((i + 1 for i, (k, _) in enumerate(quals) if k == row[0]), default=len(quals))
        quals.insert(at, row)
    return blocks


def _apply_common(part: AnnFile, doc: PlacedDoc) -> None:
    others = tuple(t for t in TARGETS if t != "ST_COMMENT")      # ST_COMMENT is the last of TARGETS
    _apply_common_targets(part, common_blocks(doc.common, others) + _st_comment_blocks(doc.common))

    if doc.source_layout == "common":
        keep = _keep_one(part, "source")
        if keep is None:
            at = len(part.rows)
            part.insert_feature(at, "COMMON" if at == 0 else "", "source", "1..E", doc.common_source)
        elif not _same_rows(part, keep, "1..E", doc.common_source):
            part.replace_feature(keep, "source", "1..E", doc.common_source)
    else:
        _delete(part, _of(part, "source"))


def _is_circular(part: AnnFile, f: AnnFeature) -> bool:
    return any("circular" in (q, v) for q, v in f.qualifiers(part.rows))


def _apply_entry(part: AnnFile, e: PlacedEntry, doc: PlacedDoc, length: int | None) -> None:
    def marker(at: int) -> str:
        return e.seq_id if at == 0 else ""

    # DBLINK: keep / rewrite the first one in place, insert below when there was none
    dblink = ([q for _, _, quals in common_blocks({"DBLINK": e.dblink}, ("DBLINK",)) for q in quals]
              if doc.dblink_layout == "entry" else [])
    new_dblink = False
    if dblink:
        keep = _keep_one(part, "DBLINK")
        if keep is None:
            new_dblink = True
        elif not _same_values(part, keep, keep.location, dblink):
            part.replace_feature(keep, "DBLINK", keep.location, dblink)
    else:
        _delete(part, _of(part, "DBLINK"))

    # TOPOLOGY: keep a circular one, insert below when there was none; "" means "not circular", so
    # only circular rows go (a TOPOLOGY linear row, which extract reads as "", stays)
    new_topology = False
    topos = _of(part, "TOPOLOGY")
    if e.topology == "circular":
        if topos and _is_circular(part, topos[0]):
            _delete(part, topos[1:])
        else:
            _delete(part, topos)
            new_topology = True
    else:
        _delete(part, [f for f in topos if _is_circular(part, f)])

    # source (layout entry without rows: the entry has no source)
    if doc.source_layout == "entry" and e.source:
        keep = _keep_one(part, "source")
        if keep is None:
            at = 0
            for f in part.features():           # after the entry's leading TOPOLOGY / DBLINK
                if f.key not in ("TOPOLOGY", "DBLINK"):
                    break
                at = f.last + 1
            part.insert_feature(at, marker(at), "source", f"1..{length}", e.source)
        else:
            location = f"1..{length}" if length is not None else keep.location
            if not _same_rows(part, keep, location, e.source):
                part.replace_feature(keep, "source", location, e.source)
    else:
        _delete(part, _of(part, "source"))

    def anchor(after_topology: bool) -> int:
        feats = part.features()
        src = next((f for f in feats if f.key == "source"), None)
        if src is not None:
            return src.first
        at = 0
        if after_topology:
            for f in feats:
                if f.key != "TOPOLOGY":
                    break
                at = f.last + 1
        return at

    if new_topology:
        at = anchor(False)
        part.insert_feature(at, marker(at), "TOPOLOGY", "", [("circular", "")])
    if new_dblink:
        at = anchor(True)
        part.insert_feature(at, marker(at), "DBLINK", "", dblink)


def apply_placed(ann: AnnFile, doc: PlacedDoc, lengths: dict[str, int]) -> None:
    """Rewrite *ann* in place from *doc* (checked by :func:`placed_issues`).

    *lengths* gives the source location ``1..<length>`` of each entry; its key order (FASTA
    order) places entries that have no rows in the .ann.
    """
    eol = ann.default_eol()
    final_newline = not ann.rows or ann.rows[-1].eol != ""
    head: list[AnnRow] = []
    parts: list[tuple[str, _Part]] = []
    for r in ann.rows:
        if r.cols[0]:
            parts.append((r.cols[0], _Part([r], eol)))
        elif parts:
            parts[-1][1].rows.append(r)
        else:
            head.append(r)
    by_id = dict(parts)
    if "COMMON" not in by_id:
        by_id["COMMON"] = _Part([], eol)
        parts.insert(0, ("COMMON", by_id["COMMON"]))
    _apply_common(by_id["COMMON"], doc)

    new: dict[str, _Part] = {}
    for e in doc.entries:
        part = by_id.get(e.seq_id)
        if part is None:
            part = new[e.seq_id] = _Part([], eol)
        _apply_entry(part, e, doc, lengths.get(e.seq_id))

    # a new entry goes before the first entry after it in FASTA order that is in the .ann
    order = list(lengths) + [e.seq_id for e in doc.entries if e.seq_id not in lengths]
    before: dict[str | None, list[str]] = {}
    nxt: str | None = None
    for sid in reversed(order):
        if sid in by_id:
            nxt = sid
        elif sid in new:
            before.setdefault(nxt, []).insert(0, sid)
    rows = list(head)
    for pid, part in parts:
        for sid in before.get(pid, []):
            rows += new[sid].rows
        rows += part.rows
    for sid in before.get(None, []):
        rows += new[sid].rows

    for r in rows[:-1]:
        if not r.eol:
            r.eol = eol
    if rows:
        rows[-1].eol = (rows[-1].eol or eol) if final_newline else ""
    ann.rows[:] = rows
