"""rRNA helpers: product and gene names of ribosomal RNA features.


Product forms accepted:
  "16S ribosomal RNA", "16S rRNA", "4.5S ribosomal RNA"
                       (size in Svedberg units + "S")
  "small subunit ribosomal RNA", "large subunit ribosomal RNA"
                       (used for organelle genes such as fungal mtDNA)

Gene names whose size can be checked against the product:
  rrn16, rrn23, rrn4.5, rrn5, ... -> 16S, 23S, 4.5S, 5S, ...
  rrs -> 16S, rrl -> 23S, rrf -> 5S   (bacterial names)
Other names (e.g. rns, rnl) are not checked.
"""

from __future__ import annotations

import re
from typing import Optional

_SIZE_PRODUCT_RE = re.compile(r"^(\d+(?:\.\d+)?)S (?:ribosomal RNA|rRNA)$")
_SUBUNIT_PRODUCT_RE = re.compile(r"^(small|large) subunit ribosomal RNA$")
_RRN_GENE_RE = re.compile(r"^rrn(\d+(?:\.\d+)?)$")
_BACTERIAL_GENES = {"rrs": "16S", "rrl": "23S", "rrf": "5S"}


def parse_rrna_product(product: str) -> Optional[str]:
    """Return a normalized product class, or None if not a known form.

    "16S ribosomal RNA" / "16S rRNA" -> "16S";
    "small subunit ribosomal RNA" -> "SSU"; "large subunit ..." -> "LSU".
    """
    m = _SIZE_PRODUCT_RE.match(product)
    if m:
        return f"{m.group(1)}S"
    m = _SUBUNIT_PRODUCT_RE.match(product)
    if m:
        return "SSU" if m.group(1) == "small" else "LSU"
    return None


def rrna_gene_size(gene: str) -> Optional[str]:
    """Return the size ("16S", "4.5S", ...) implied by a gene name, or None."""
    m = _RRN_GENE_RE.match(gene)
    if m:
        return f"{m.group(1)}S"
    return _BACTERIAL_GENES.get(gene)
