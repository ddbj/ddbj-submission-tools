"""tRNA helpers: amino acid names, gene names and the anticodon qualifier.


Conventions used in the samples:
  gene      trnH-GUG   trn + one-letter amino acid (fM for initiator Met)
                       + "-" + anticodon written in RNA (5'->3')
  product   tRNA-His   "tRNA-" + three-letter amino acid
  anticodon (pos:106..108,aa:Met)   pos is a location on the entry
                       sequence, possibly complement(...)
"""

from __future__ import annotations

import re
from typing import Optional, Tuple

# One-letter -> three-letter amino acid codes used in tRNA names
AA1_TO_AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "U": "Sec", "O": "Pyl",
    "fM": "Met",                     # initiator tRNA (trnfM)
}

_GENE_RE = re.compile(r"^trn(fM|[A-Z])(?:-([ACGU]{3}))?(?:\(.*\))?$")
_GENE_PAREN_RE = re.compile(r"^trn(fM|[A-Z])\(([acgtu]{3})\)$")
_PRODUCT_RE = re.compile(r"^tRNA-([A-Z][a-z]{2})$")
_ANTICODON_RE = re.compile(r"^\(pos:(.+),aa:([A-Za-z]+)(?:,seq:([acgtu]{3}))?\)$")


def parse_trna_gene(gene: str) -> Optional[Tuple[str, Optional[str]]]:
    """Parse a tRNA gene name.

    Returns (three-letter amino acid, anticodon in DNA upper case or None),
    or None if the name does not follow the trnX-NNN / trnX(nnn) style.
    Examples: "trnH-GUG" -> ("His", "GTG"); "trnM(cat)" -> ("Met", "CAT");
    "trnfM-CAU" -> ("Met", "CAT"); "trnL" -> ("Leu", None).
    """
    m = _GENE_PAREN_RE.match(gene)
    if m:
        return AA1_TO_AA3[m.group(1)], m.group(2).upper().replace("U", "T")
    m = _GENE_RE.match(gene)
    if not m or m.group(1) not in AA1_TO_AA3:
        return None
    anticodon = m.group(2).replace("U", "T") if m.group(2) else None
    return AA1_TO_AA3[m.group(1)], anticodon


def parse_trna_product(product: str) -> Optional[str]:
    """Return the three-letter amino acid of "tRNA-Xxx", or None."""
    m = _PRODUCT_RE.match(product)
    return m.group(1) if m else None


def parse_anticodon(value: str) -> Optional[Tuple[str, str]]:
    """Parse the anticodon qualifier "(pos:LOCATION,aa:Xxx[,seq:nnn])".

    Returns (location text, amino acid), or None if malformed.
    """
    m = _ANTICODON_RE.match(value)
    if not m:
        return None
    return m.group(1), m.group(2)
