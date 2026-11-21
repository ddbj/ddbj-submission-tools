"""Reader for DDBJ MSS sequence files (.fasta).


File structure (summary):
  - Each entry: a header line ">entry_name", sequence lines, and a
    terminator line "//".
  - Entry names correspond to the Entry column of the .ann file
    (same names, same order).

Sequence files can be larger than 1 GB, so this module reads them line by
line. parse() keeps only the entry name and sequence length;
read_sequence() / iter_sequences() return the sequence itself.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterator, List, Optional, Union

PathLike = Union[str, Path]

TERMINATOR = "//"


class FastaParseError(Exception):
    """Raised when the sequence file violates the structural rules."""

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


@dataclass
class SeqRecordInfo:
    """Summary of one sequence entry (the sequence itself is not kept)."""
    name: str
    length: int
    line_no: int          # line number of the ">" header

    def __repr__(self) -> str:
        return f"<SeqRecordInfo {self.name} length={self.length}>"


def parse(path: PathLike) -> Iterator[SeqRecordInfo]:
    """Iterate over entries, yielding name and length for each.

    Raises FastaParseError for sequence lines outside an entry, a header
    inside an unterminated entry, or a missing final "//".
    """
    name: Optional[str] = None
    header_line = 0
    length = 0

    with open(path, encoding="ascii", errors="replace") as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    raise FastaParseError(
                        f"entry '{name}' is not terminated by '//'",
                        path, line_no)
                name = line[1:].strip()
                if not name:
                    raise FastaParseError("empty entry name", path, line_no)
                header_line = line_no
                length = 0
            elif line.strip() == TERMINATOR:
                if name is None:
                    raise FastaParseError(
                        "'//' without a preceding entry", path, line_no)
                yield SeqRecordInfo(name=name, length=length,
                                    line_no=header_line)
                name = None
            elif line.strip() == "":
                # Blank lines are not allowed inside an entry
                if name is not None:
                    raise FastaParseError(
                        f"blank line inside entry '{name}'", path, line_no)
            else:
                if name is None:
                    raise FastaParseError(
                        "sequence line outside an entry", path, line_no)
                length += len(line.strip())

    if name is not None:
        raise FastaParseError(
            f"entry '{name}' is not terminated by '//' (end of file)",
            path, header_line)


def read_names(path: PathLike) -> List[str]:
    """Return entry names in file order."""
    return [rec.name for rec in parse(path)]


def find(path: PathLike, name: str) -> Optional[SeqRecordInfo]:
    """Return the record for the given entry name, or None.

    Reading stops as soon as the entry is found.
    """
    for rec in parse(path):
        if rec.name == name:
            return rec
    return None


def iter_sequences(path: PathLike) -> Iterator[tuple]:
    """Iterate over (name, sequence) pairs in file order.

    Each sequence is held in memory only while it is being yielded.
    Structural errors are reported as in parse().
    """
    name: Optional[str] = None
    header_line = 0
    chunks: List[str] = []

    with open(path, encoding="ascii", errors="replace") as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    raise FastaParseError(
                        f"entry '{name}' is not terminated by '//'",
                        path, line_no)
                name = line[1:].strip()
                if not name:
                    raise FastaParseError("empty entry name", path, line_no)
                header_line = line_no
                chunks = []
            elif line.strip() == TERMINATOR:
                if name is None:
                    raise FastaParseError(
                        "'//' without a preceding entry", path, line_no)
                yield name, "".join(chunks)
                name = None
                chunks = []
            elif line.strip() == "":
                if name is not None:
                    raise FastaParseError(
                        f"blank line inside entry '{name}'", path, line_no)
            else:
                if name is None:
                    raise FastaParseError(
                        "sequence line outside an entry", path, line_no)
                chunks.append(line.strip())

    if name is not None:
        raise FastaParseError(
            f"entry '{name}' is not terminated by '//' (end of file)",
            path, header_line)


def read_sequence(path: PathLike, name: str) -> Optional[str]:
    """Return the sequence of the given entry, or None if not found.

    Reading stops as soon as the entry is read.
    """
    for entry_name, seq in iter_sequences(path):
        if entry_name == name:
            return seq
    return None


def rename_entries(src: PathLike, dst: PathLike,
                   rename: "Callable[[str], str]") -> int:
    """Copy a sequence file, renaming entries with rename(name) -> new name.

    Only ">" lines change; sequence lines are copied as they are, line by
    line, so large files are handled. Returns the number of renamed entries.
    """
    n = 0
    with open(src, encoding="ascii", errors="replace", newline="") as fin, \
            open(dst, "w", encoding="ascii", newline="") as fout:
        for raw in fin:
            if raw.startswith(">"):
                name = raw[1:].strip()
                new = rename(name)
                if new != name:
                    n += 1
                eol = "\r\n" if raw.endswith("\r\n") else "\n"
                fout.write(f">{new}{eol}")
            else:
                fout.write(raw)
    return n
