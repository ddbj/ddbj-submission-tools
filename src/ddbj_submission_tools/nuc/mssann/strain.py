"""Helpers for the strain qualifier: culture collection numbers.


When a strain value starts with a culture collection code, the code and
the number must be separated by a space ("JCM11017" -> "JCM 11017").
Culture collection codes come from NCBI's coll_dump.txt
(https://ftp.ncbi.nih.gov/pub/taxonomy/coll_dump.txt): column 1 is the
code, column 2 the collection type letters; codes whose column 2 contains
"c" (culture collection) are used. Codes in the separator exception list
(mssann/data/strain_separator_exceptions.tsv) use another separator,
e.g. NIES-1234.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Union

PathLike = Union[str, Path]

DEFAULT_EXCEPTIONS = (Path(__file__).resolve().parent / "data"
                      / "strain_separator_exceptions.tsv")
DEFAULT_SEPARATOR = " "


class StrainDataError(ValueError):
    """Raised when coll_dump.txt or the exception list is unusable."""


def load_culture_collection_codes(path: PathLike) -> List[str]:
    """Return culture collection codes (column 2 contains "c"), longest first.

    Longest first, so that a longer code wins over a shorter one that is
    its prefix (e.g. "CCUG" before "CC").
    """
    path = Path(path)
    codes = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line_no, raw in enumerate(fh, start=1):
            cols = raw.rstrip("\r\n").split("\t")
            if len(cols) < 2 or not cols[0].strip():
                continue
            if "c" in cols[1]:
                codes.add(cols[0].strip())
    if not codes:
        raise StrainDataError(f"{path}: no culture collection codes found "
                              "(expected column 2 containing 'c')")
    return sorted(codes, key=lambda c: (-len(c), c))


def load_separator_exceptions(path: Optional[PathLike] = None) -> Dict[str, str]:
    """Load {code: separator} from the exception list."""
    path = Path(path) if path is not None else DEFAULT_EXCEPTIONS
    result: Dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\r\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2 or not cols[0] or not cols[1]:
                raise StrainDataError(f"{path}:{line_no}: need code and "
                                      "separator columns")
            if cols[0] in result:
                raise StrainDataError(f"{path}:{line_no}: duplicated code "
                                      f"'{cols[0]}'")
            result[cols[0]] = cols[1]
    return result


@dataclass
class StrainCheck:
    """Result of checking one strain value."""
    value: str
    ok: bool
    code: Optional[str] = None       # culture collection code found
    message: str = ""
    fix: Optional[str] = None


def check_strain_value(value: str, codes: List[str],
                       exceptions: Dict[str, str]) -> StrainCheck:
    """Check the separator between a culture collection code and number.

    The value is checked only when it starts with a code (case-sensitive)
    that is directly followed by a digit or by a separator and a digit.
    The expected separator is a space, or the one in the exception list.
      code + digit                -> error, separator missing
      code + other sep + digit    -> error for exception codes (NIES 12)
                                     and space-codes written with the
                                     exception style are left unchecked
    """
    for code in codes:
        if not value.startswith(code):
            continue
        rest = value[len(code):]
        if not rest:
            continue
        expected = exceptions.get(code, DEFAULT_SEPARATOR)
        if rest[0].isdigit():
            fix = code + expected + rest
            return StrainCheck(value, False, code,
                               f"separator '{expected}' missing after "
                               f"culture collection code {code}", fix)
        if rest.startswith(expected) and rest[len(expected):][:1].isdigit():
            return StrainCheck(value, True, code, f"{code}{expected}number")
        if code in exceptions and rest[:1] in (" ", "-", "_") \
                and rest[1:2].isdigit():
            fix = code + expected + rest[1:]
            return StrainCheck(value, False, code,
                               f"{code} uses '{expected}' as separator", fix)
        # Followed by a letter etc.: this code does not apply; try others
    return StrainCheck(value, True, None, "no culture collection code")
