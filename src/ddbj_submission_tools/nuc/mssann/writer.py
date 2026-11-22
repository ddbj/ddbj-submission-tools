"""Write parsed entries back to the MSS annotation format (.ann).


Each qualifier becomes one line of 5 tab-separated columns. The Entry
column is written only on the first line of an entry, Feature and
Location only on the first line of a feature (the inverse of the parser).
A feature with no qualifier at all (5'UTR, 3'UTR, D-loop, centromere,
telomere and other region-only features) becomes a single line whose
Qualifier and Value columns are empty, as long as it has a location.
Lines end with LF. Reading a file and writing it unchanged reproduces it
exactly, except that CRLF becomes LF and spreadsheet quoting "..." is
removed (use keep_quotes=True to write quoted fields back as they were).
"""

from __future__ import annotations

from typing import IO, Iterable, List, Tuple

from .annotation import Entry, Feature


class WriteError(ValueError):
    """Raised when an entry cannot be written in the .ann format."""


def _check(text: str, what: str, line_no: int) -> str:
    if "\t" in text or "\n" in text or "\r" in text:
        raise WriteError(f"line {line_no}: {what} contains a tab or newline: "
                         f"{text!r}")
    return text


def _quote(text: str) -> str:
    return '"' + text.replace('"', '""') + '"'


def format_entry(entry: Entry, keep_quotes: bool = False
                 ) -> Tuple[List[str], List[Feature]]:
    """Return (lines without newline, features dropped).

    A feature with no qualifier is written as one line with empty
    Qualifier and Value columns when it has a location. Without a
    location such a feature carries nothing at all, so it is dropped and
    returned for the caller to report.
    """
    lines: List[str] = []
    dropped: List[Feature] = []
    first_line_of_entry = True
    for feat in entry.features:
        if not feat.qualifiers and not feat.location:
            dropped.append(feat)
            continue
        location = feat.location or ""
        if keep_quotes and feat.location_quoted:
            location = _quote(location)

        if not feat.qualifiers:
            cols = [
                entry.name if first_line_of_entry else "",
                feat.key,
                location,
                "",
                "",
            ]
            for name, text in zip(("entry name", "feature key", "location"),
                                  cols):
                _check(text, name, feat.line_no)
            lines.append("\t".join(cols))
            first_line_of_entry = False
            continue

        for i, q in enumerate(feat.qualifiers):
            value = _quote(q.value) if keep_quotes and q.quoted else q.value
            cols = [
                entry.name if first_line_of_entry else "",
                feat.key if i == 0 else "",
                location if i == 0 else "",
                q.key,
                value,
            ]
            for name, text in zip(("entry name", "feature key", "location",
                                   "qualifier", "value"), cols):
                _check(text, name, q.line_no)
            lines.append("\t".join(cols))
            first_line_of_entry = False
    return lines, dropped


def write_entries(entries: Iterable[Entry], fh: IO[str],
                  keep_quotes: bool = False) -> List[Tuple[Entry, Feature]]:
    """Write entries to an open text file; return dropped (entry, feature)."""
    dropped_all: List[Tuple[Entry, Feature]] = []
    for entry in entries:
        lines, dropped = format_entry(entry, keep_quotes=keep_quotes)
        if not lines:
            raise WriteError(f"entry '{entry.name}' has no line left; "
                             "an entry needs at least one line")
        for line in lines:
            fh.write(line + "\n")
        dropped_all.extend((entry, f) for f in dropped)
    return dropped_all
