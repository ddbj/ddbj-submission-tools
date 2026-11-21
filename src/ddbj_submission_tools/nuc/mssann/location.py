"""Parser for feature locations in DDBJ MSS annotation files (.ann).


Supported syntax (1-based, both ends inclusive):
    n                       single base (start == end is always written as "n")
    a^b                     site between bases a and b (length 0); b == a + 1,
                            or "seq_length^1" on a circular sequence
    a..b                    range; a < b is required
    <a..b  a..>b  <a..>b    partial start / partial end / both
    complement(X)           X on the complementary strand
    join(X,Y,...)           X, Y, ... joined into one contiguous feature
    order(X,Y,...)          X, Y, ... in that order, not necessarily joined
    join(complement(X),Y)   mixed-strand join (e.g. trans-splicing)
    E                       meta position "last base" (COMMON entry only);
                            resolved when seq_length is given

Nested forms such as complement(join(...)) and join(complement(...),...)
are accepted. The "a.b" (one base within a range) form is not supported.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import List, Optional

META_END = "E"

# One position: optional partial marker, then digits or the meta "E".
# Two positions are separated by ".." (range) or "^" (between bases).
_POS_RE = r"([<>]?)(\d+|E)"
_SPAN_RE = re.compile(rf"^{_POS_RE}(?:(\.\.|\^){_POS_RE})?$")
_SPAN_TOKEN_RE = re.compile(r"[<>]?(?:\d+|E)(?:(?:\.\.|\^)[<>]?(?:\d+|E))?")


class LocationError(ValueError):
    """Raised when a location string cannot be parsed or is invalid."""


# Complement table for IUPAC nucleotide codes (case preserved)
_COMPLEMENT = str.maketrans(
    "ACGTUMRWSYKVHDBNacgtumrwsykvhdbn",
    "TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn",
)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a nucleotide sequence (IUPAC)."""
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class Span:
    """One contiguous region."""
    start: int
    end: int
    partial_start: bool = False   # "<" before start
    partial_end: bool = False     # ">" before end
    complement: bool = False      # on the complementary strand
    between: bool = False         # "a^b": site between two bases

    @property
    def length(self) -> int:
        """Number of bases; 0 for a site between bases."""
        if self.between:
            return 0
        return self.end - self.start + 1

    @property
    def wraps_origin(self) -> bool:
        """True for "seq_length^1" on a circular sequence."""
        return self.between and self.end < self.start

    def __str__(self) -> str:
        if self.between:
            s = f"{self.start}^{self.end}"
        elif self.start == self.end:
            mark = "<" if self.partial_start else (">" if self.partial_end else "")
            s = f"{mark}{self.start}"
        else:
            s = (("<" if self.partial_start else "") + str(self.start) + ".."
                 + (">" if self.partial_end else "") + str(self.end))
        return f"complement({s})" if self.complement else s


@dataclass
class Location:
    """A parsed location string."""
    raw: str
    spans: List[Span] = field(default_factory=list)
    operator: Optional[str] = None   # "join", "order" or None
    complement: bool = False         # outermost complement(...)

    @property
    def start(self) -> int:
        """Smallest coordinate over all spans."""
        return min(s.start for s in self.spans)

    @property
    def end(self) -> int:
        """Largest coordinate over all spans."""
        return max(s.end for s in self.spans)

    @property
    def length(self) -> int:
        """Total number of bases covered by the spans (0 for "a^b")."""
        return sum(s.length for s in self.spans)

    @property
    def strand(self) -> int:
        """+1 if all spans are forward, -1 if all complement, 0 if mixed."""
        flags = {s.complement for s in self.spans}
        if flags == {False}:
            return 1
        if flags == {True}:
            return -1
        return 0

    @property
    def is_between(self) -> bool:
        """True if the location is a single site between two bases."""
        return len(self.spans) == 1 and self.spans[0].between

    @property
    def is_partial(self) -> bool:
        return any(s.partial_start or s.partial_end for s in self.spans)

    def product_spans(self) -> List[Span]:
        """Return spans in product order (5' -> 3' of the feature).

        This is the written order, except for an outer complement(join(...))
        whose spans are taken in reverse order.
        """
        if self.complement and self.operator is not None:
            return list(reversed(self.spans))
        return list(self.spans)

    @property
    def five_prime_partial(self) -> bool:
        """True if the 5' end of the product is partial.

        On the forward strand this is '<' of the first span; on the
        complementary strand it is '>' of that span (the higher coordinate).
        """
        first = self.product_spans()[0]
        return first.partial_end if first.complement else first.partial_start

    @property
    def three_prime_partial(self) -> bool:
        """True if the 3' end of the product is partial."""
        last = self.product_spans()[-1]
        return last.partial_start if last.complement else last.partial_end

    def extract(self, seq: str) -> str:
        """Return the nucleotide sequence described by this location.

        seq is the whole entry sequence (position 1 = seq[0]).
        Spans are concatenated in the written order (5' -> 3' of the
        product). For an outer complement(join(...)) the joined sequence is
        reverse-complemented, i.e. the spans are taken in reverse order.
        A site between bases ("a^b") yields an empty string.
        """
        # complement(join(A,B)) == revcomp(A+B) == revcomp(B)+revcomp(A)
        parts = []
        for span in self.product_spans():
            if span.between:
                continue
            if span.end > len(seq):
                raise LocationError(
                    f"span {span} exceeds sequence length {len(seq)}")
            sub = seq[span.start - 1:span.end]
            parts.append(reverse_complement(sub) if span.complement else sub)
        return "".join(parts)

    def __str__(self) -> str:
        return self.raw


# ---------------------------------------------------------------------------
# Recursive-descent parser
# ---------------------------------------------------------------------------

class _Parser:
    def __init__(self, text: str, seq_length: Optional[int], circular: bool):
        self.text = text
        self.pos = 0
        self.seq_length = seq_length
        self.circular = circular

    def error(self, message: str) -> LocationError:
        return LocationError(f"{message} at position {self.pos + 1} "
                             f"in '{self.text}'")

    def peek_word(self, word: str) -> bool:
        return self.text.startswith(word + "(", self.pos)

    def expect(self, char: str) -> None:
        if self.pos >= len(self.text) or self.text[self.pos] != char:
            raise self.error(f"expected '{char}'")
        self.pos += 1

    def parse_element(self, complement: bool) -> tuple:
        """Parse one element. Return (spans, operator) of that element."""
        if self.peek_word("complement"):
            self.pos += len("complement(")
            spans, op = self.parse_element(not complement)
            self.expect(")")
            return spans, op
        for op in ("join", "order"):
            if self.peek_word(op):
                self.pos += len(op) + 1
                spans: List[Span] = []
                while True:
                    sub, sub_op = self.parse_element(complement)
                    if sub_op is not None:
                        raise self.error(f"nested {sub_op}() inside {op}()")
                    if any(sp.between for sp in sub):
                        raise self.error(f"'a^b' cannot be used inside {op}()")
                    spans.extend(sub)
                    if self.pos < len(self.text) and self.text[self.pos] == ",":
                        self.pos += 1
                        continue
                    break
                self.expect(")")
                return spans, op
        return [self.parse_span(complement)], None

    def parse_span(self, complement: bool) -> Span:
        """Parse a simple span such as '<1..>200', '1010' or '123^124'."""
        m = _SPAN_TOKEN_RE.match(self.text, self.pos)
        if not m:
            raise self.error("expected a position")
        token = m.group(0)
        self.pos += len(token)
        sm = _SPAN_RE.match(token)
        mark1, pos1, sep, mark2, pos2 = sm.groups()

        if sep == "^":
            return self.parse_between(token, mark1, pos1, mark2, pos2,
                                      complement)

        start = self.resolve(pos1)
        if pos2 is None:
            # Single base: "n", "<n" or ">n"
            return Span(start, start, partial_start=(mark1 == "<"),
                        partial_end=(mark1 == ">"), complement=complement)

        end = self.resolve(pos2)
        if mark1 == ">":
            raise self.error("'>' is not allowed before the start position")
        if mark2 == "<":
            raise self.error("'<' is not allowed before the end position")
        if start == end:
            raise self.error(f"'{token}' has start == end; write '{start}'")
        if start > end:
            raise self.error(f"'{token}' has start > end")
        return Span(start, end, partial_start=(mark1 == "<"),
                    partial_end=(mark2 == ">"), complement=complement)

    def parse_between(self, token: str, mark1: str, pos1: str, mark2: str,
                      pos2: str, complement: bool) -> Span:
        """Parse 'a^b' (site between bases a and b)."""
        if mark1 or mark2:
            raise self.error(f"'{token}': '<' or '>' cannot be used with '^'")
        a = self.resolve(pos1)
        b = self.resolve(pos2)
        if b == a + 1:
            return Span(a, b, complement=complement, between=True)
        if b == 1:
            # Only "last base ^ first base" on a circular sequence
            if not self.circular:
                raise self.error(
                    f"'{token}': 'n^1' is allowed only on a circular sequence")
            if self.seq_length is None:
                raise self.error(
                    f"'{token}': 'n^1' needs the sequence length to check n")
            if a != self.seq_length:
                raise self.error(
                    f"'{token}': on a circular sequence only "
                    f"'{self.seq_length}^1' may wrap the origin")
            return Span(a, b, complement=complement, between=True)
        raise self.error(f"'{token}': '^' requires b == a + 1")

    def resolve(self, value: str) -> int:
        """Convert a position string to int, resolving the meta 'E'."""
        if value == META_END:
            if self.seq_length is None:
                raise self.error("meta position 'E' needs the sequence length")
            n = self.seq_length
        else:
            n = int(value)
        if n < 1:
            raise self.error("position must be >= 1")
        if self.seq_length is not None and n > self.seq_length:
            raise self.error(
                f"position {n} exceeds sequence length {self.seq_length}")
        return n


def parse_location(text: str, seq_length: Optional[int] = None,
                   circular: bool = False) -> Location:
    """Parse a location string.

    seq_length (optional) resolves the meta position "E" and enables the
    check that every position is within the sequence.
    circular=True allows "seq_length^1" (site between the last and the
    first base); use mssann.is_circular() to get it from TOPOLOGY.
    """
    if text is None or text == "":
        raise LocationError("empty location")
    if any(c.isspace() for c in text):
        raise LocationError(f"whitespace is not allowed in '{text}'")

    parser = _Parser(text, seq_length, circular)
    outer_complement = parser.peek_word("complement")
    spans, op = parser.parse_element(False)
    if parser.pos != len(text):
        raise parser.error("unexpected trailing text")

    loc = Location(raw=text, spans=spans, operator=op,
                   complement=outer_complement)
    return loc
