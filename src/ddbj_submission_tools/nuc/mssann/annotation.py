"""Parser and data model for DDBJ MSS annotation files (.ann).


File structure (summary):
  - Each line has 5 tab-separated columns:
        Entry, Feature, Location, Qualifier, Value
  - The Entry column is filled only on the first line of an entry.
  - The Feature column is filled only on the first line of a feature.
  - The Location column is filled only on the first line of a feature
    (biological features and PRIMARY_CONTIG).
  - One Qualifier/Value pair per line; Value may be empty
    (e.g. environmental_sample, TOPOLOGY/circular). A feature may also
    carry no qualifier at all, in which case its single line has empty
    Qualifier and Value columns (5'UTR, 3'UTR, D-loop, centromere,
    telomere ...).
  - A blank line marks the end of the annotation data.
  - Spreadsheet software may save a field as "..." (CSV-style quoting,
    "" for a literal quote). The parser removes such quoting by default
    and records it in Feature.location_quoted / Qualifier.quoted.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, List, Optional, Union

PathLike = Union[str, Path]

NUM_COLUMNS = 5
COMMON = "COMMON"

# DDBJ-specific features (no Location except PRIMARY_CONTIG)
DDBJ_FEATURES = frozenset({
    "SUBMITTER", "REFERENCE", "DATE", "COMMENT", "ST_COMMENT", "DIVISION",
    "DATATYPE", "KEYWORD", "DBLINK", "TOPOLOGY", "PRIMARY_CONTIG",
})


# Allowed values of the pseudogene qualifier. The old "pseudo" qualifier
# cannot be used in new submissions.
PSEUDOGENE_VALUES = ("processed", "unprocessed", "unitary", "allelic",
                     "unknown")


class AnnParseError(Exception):
    """Raised when the file violates the structural rules of the spec."""

    def __init__(self, message: str, path: Optional[PathLike] = None,
                 line_no: Optional[int] = None):
        self.path = str(path) if path is not None else None
        self.line_no = line_no
        where = ""
        if self.path is not None:
            where += self.path
        if line_no is not None:
            where += f":{line_no}"
        super().__init__(f"{where}: {message}" if where else message)


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class Qualifier:
    """One Qualifier/Value line."""
    key: str
    value: str
    line_no: int
    quoted: bool = False      # value was written as "..." in the file


@dataclass
class Feature:
    """A feature (biological or DDBJ-specific) with its qualifiers."""
    key: str
    location: Optional[str]
    line_no: int
    qualifiers: List[Qualifier] = field(default_factory=list)
    location_quoted: bool = False   # location was written as "..."

    @property
    def is_ddbj_feature(self) -> bool:
        return self.key in DDBJ_FEATURES

    def get_all(self, key: str) -> List[str]:
        """Return all values of a qualifier in file order (duplicates kept)."""
        return [q.value for q in self.qualifiers if q.key == key]

    def get(self, key: str, default: Optional[str] = None) -> Optional[str]:
        """Return the first value of a qualifier, or default if absent."""
        for q in self.qualifiers:
            if q.key == key:
                return q.value
        return default

    def has(self, key: str) -> bool:
        """Return True if the qualifier exists (even with an empty value)."""
        return any(q.key == key for q in self.qualifiers)

    @property
    def is_pseudo(self) -> bool:
        """True if the feature is marked as a pseudogene (pseudogene or the
        old pseudo qualifier)."""
        return self.has("pseudogene") or self.has("pseudo")

    def pseudo_problems(self) -> List[str]:
        """Return messages for pseudo / pseudogene qualifiers that do not
        follow the current submission rules (empty list if none)."""
        problems = []
        allowed = ", ".join(PSEUDOGENE_VALUES)
        if self.has("pseudo"):
            problems.append("'pseudo' cannot be used in new submissions; "
                            f"use 'pseudogene' ({allowed})")
        values = self.get_all("pseudogene")
        if len(values) > 1:
            problems.append(f"pseudogene is given {len(values)} times")
        for value in values:
            if value not in PSEUDOGENE_VALUES:
                problems.append(f"invalid pseudogene value '{value}' "
                                f"(allowed: {allowed})")
        return problems

    def __repr__(self) -> str:
        loc = f" {self.location}" if self.location else ""
        return f"<Feature {self.key}{loc} ({len(self.qualifiers)} qualifiers)>"


@dataclass
class Entry:
    """An entry: COMMON or a sequence entry."""
    name: str
    line_no: int
    features: List[Feature] = field(default_factory=list)

    @property
    def is_common(self) -> bool:
        return self.name == COMMON

    def get_features(self, key: str) -> List[Feature]:
        """Return all features with the given key in file order."""
        return [f for f in self.features if f.key == key]

    def get_feature(self, key: str) -> Optional[Feature]:
        """Return the first feature with the given key, or None."""
        for f in self.features:
            if f.key == key:
                return f
        return None

    def get_values(self, feature_key: str, qualifier_key: str) -> List[str]:
        """Return all values of qualifier_key across all feature_key features."""
        values: List[str] = []
        for f in self.get_features(feature_key):
            values.extend(f.get_all(qualifier_key))
        return values

    def __repr__(self) -> str:
        return f"<Entry {self.name} ({len(self.features)} features)>"


@dataclass
class AnnotationFile:
    """Whole annotation file: optional COMMON entry plus sequence entries."""
    path: Optional[str]
    common: Optional[Entry]
    entries: List[Entry]

    def __iter__(self) -> Iterator[Entry]:
        return iter(self.entries)

    def __len__(self) -> int:
        return len(self.entries)

    def get_entry(self, name: str) -> Optional[Entry]:
        for e in self.entries:
            if e.name == name:
                return e
        return None

    def entry_names(self) -> List[str]:
        return [e.name for e in self.entries]


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def split_row(line: str, path: Optional[PathLike] = None,
              line_no: Optional[int] = None) -> List[str]:
    """Split one line into exactly 5 columns.

    Missing trailing columns are padded with "" (some editors strip trailing
    tabs). Extra columns are dropped if they are empty (trailing tabs);
    extra columns with content are a structural error.
    """
    cols = line.rstrip("\r\n").split("\t")
    if len(cols) > NUM_COLUMNS:
        # Extra empty columns (trailing tabs left by spreadsheet software)
        # are dropped; extra columns with content are an error.
        if any(c.strip() for c in cols[NUM_COLUMNS:]):
            raise AnnParseError(
                f"too many columns ({len(cols)}; expected {NUM_COLUMNS}); "
                "columns 6 and beyond are not empty", path, line_no)
        cols = cols[:NUM_COLUMNS]
    if len(cols) < NUM_COLUMNS:
        cols += [""] * (NUM_COLUMNS - len(cols))
    return cols


def unquote(text: str) -> tuple:
    """Remove CSV-style quoting added by spreadsheet software.

    A field wrapped in double quotes loses the outer quotes and each
    doubled quote becomes one quote (e.g. "a,b" -> a,b), returning
    (text, True). Anything else is returned unchanged with False.
    """
    if len(text) >= 2 and text[0] == '"' and text[-1] == '"':
        return text[1:-1].replace('""', '"'), True
    return text, False


def _iter_lines(path: PathLike, check_trailing: bool) -> Iterator[tuple]:
    """Yield (line_no, line) up to the first blank line.

    If check_trailing is True, any non-blank line after the first blank
    line raises AnnParseError (blank lines are allowed only at the end of an .ann file).
    """
    with open(path, encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            if line.strip() == "":
                if check_trailing:
                    for extra_no, extra in enumerate(fh, start=line_no + 1):
                        if extra.strip() != "":
                            raise AnnParseError(
                                "data found after a blank line "
                                f"(blank line at {line_no})", path, extra_no)
                return
            yield line_no, line


def parse(path: PathLike, check_trailing: bool = True,
          strip_quotes: bool = True) -> Iterator[Entry]:
    """Iterate over entries (COMMON first if present) one at a time.

    Memory use is bounded by the size of a single entry, so this is
    suitable for genome-scale files. Stop iterating early to avoid reading
    the rest of the file (trailing-data check is then skipped).
    strip_quotes removes spreadsheet "..." quoting from Location and Value.
    """
    entry: Optional[Entry] = None
    feature: Optional[Feature] = None
    seen_names = set()

    for line_no, line in _iter_lines(path, check_trailing):
        e_col, f_col, loc_col, q_col, v_col = split_row(line, path, line_no)
        loc_quoted = v_quoted = False
        if strip_quotes:
            loc_col, loc_quoted = unquote(loc_col)
            v_col, v_quoted = unquote(v_col)

        if e_col:
            if not f_col:
                raise AnnParseError(
                    f"entry '{e_col}' starts without a feature", path, line_no)
            if e_col in seen_names:
                raise AnnParseError(
                    f"duplicate entry name '{e_col}'", path, line_no)
            if e_col == COMMON and seen_names:
                raise AnnParseError(
                    "COMMON must be the first entry", path, line_no)
            if entry is not None:
                yield entry
            seen_names.add(e_col)
            entry = Entry(name=e_col, line_no=line_no)
            feature = None
        elif entry is None:
            raise AnnParseError(
                "first line has an empty Entry column", path, line_no)

        if f_col:
            feature = Feature(key=f_col, location=loc_col or None,
                              line_no=line_no, location_quoted=loc_quoted)
            entry.features.append(feature)
        else:
            if feature is None:
                raise AnnParseError(
                    "qualifier line before any feature", path, line_no)
            if loc_col:
                raise AnnParseError(
                    "Location given on a line without a Feature",
                    path, line_no)

        if not q_col:
            # A feature may carry no qualifier at all (5'UTR, 3'UTR,
            # D-loop, centromere, telomere ...). Only a line that
            # continues a feature must have one.
            if not f_col:
                raise AnnParseError("empty Qualifier column", path, line_no)
            if v_col:
                raise AnnParseError(
                    "Value given without a Qualifier", path, line_no)
            continue
        feature.qualifiers.append(
            Qualifier(key=q_col, value=v_col, line_no=line_no,
                      quoted=v_quoted))

    if entry is not None:
        yield entry


def read(path: PathLike, check_trailing: bool = True,
         strip_quotes: bool = True) -> AnnotationFile:
    """Read the whole file into an AnnotationFile object."""
    common: Optional[Entry] = None
    entries: List[Entry] = []
    for entry in parse(path, check_trailing=check_trailing,
                       strip_quotes=strip_quotes):
        if entry.is_common:
            common = entry
        else:
            entries.append(entry)
    return AnnotationFile(path=str(path), common=common, entries=entries)


def read_common(path: PathLike) -> Optional[Entry]:
    """Return only the COMMON entry (or None), reading as little as possible."""
    for entry in parse(path, check_trailing=False):
        return entry if entry.is_common else None
    return None


def find_entry(path: PathLike, name: str) -> tuple:
    """Return (common, entry) for the given entry name.

    common is the COMMON entry or None; entry is None if the name is not
    in the .ann file. Reading stops as soon as the entry is found.
    """
    common: Optional[Entry] = None
    for entry in parse(path, check_trailing=False):
        if entry.is_common:
            common = entry
        elif entry.name == name:
            return common, entry
    return common, None


def effective_source(entry: Optional[Entry],
                     common: Optional[Entry]) -> Optional[Feature]:
    """Return the source feature that applies to an entry.

    The entry's own source takes precedence. Otherwise the COMMON source
    (written with location "1..E") applies to every entry, including
    entries that appear only in the .fasta file (entry is None).
    """
    if entry is not None:
        src = entry.get_feature("source")
        if src is not None:
            return src
    if common is not None:
        return common.get_feature("source")
    return None


def is_circular(entry: Optional[Entry], common: Optional[Entry]) -> bool:
    """Return True if TOPOLOGY/circular applies to the entry.

    TOPOLOGY may be written in the entry itself or in COMMON.
    """
    for e in (entry, common):
        if e is None:
            continue
        for f in e.get_features("TOPOLOGY"):
            if f.has("circular"):
                return True
    return False
