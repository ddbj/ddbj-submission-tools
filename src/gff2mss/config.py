from __future__ import annotations

import tomllib
from dataclasses import dataclass, field

from ddbj_gff.errors import Diagnostic, GffParseError, Severity

_KNOWN_SECTIONS = {"locus_tag", "cds", "product", "transcript"}
# Metadata now comes from the common JSON (and sequence roles), as in mss_builder.
_MOVED_TO_COMMON = {"source": "SOURCE", "assembly_gap": "ASSEMBLY_GAP"}


@dataclass
class MssConfig:
    """Annotation conversion settings. Metadata (source, gaps) is in the common JSON."""
    locus_tag_prefix: str = ""
    locus_tag_width: int = 6
    locus_tag_start: int = 10
    locus_tag_step: int = 10
    transl_table: int = 1
    product_default: str = "hypothetical protein"
    product_map_path: str | None = None
    product_map: dict = field(default_factory=dict)
    transcript_mode: str = "nonredundant"
    emit_mrna: bool = True   # organelle/prokaryotic profiles set False (CDS only, no mRNA feature)


def load_config(path: str) -> tuple[MssConfig, list[Diagnostic]]:
    with open(path, "rb") as fh:
        data = tomllib.load(fh)
    diags: list[Diagnostic] = []
    for key in data:
        if key in _MOVED_TO_COMMON:
            diags.append(Diagnostic(Severity.WARNING, None, "moved-to-common",
                                    f"[{key}] in the gff2mss config is ignored; put it in the "
                                    f"common JSON ({_MOVED_TO_COMMON[key]})"))
        elif key not in _KNOWN_SECTIONS:
            diags.append(Diagnostic(Severity.WARNING, None, "unknown-config-key",
                                    f"unknown config section [{key}] ignored"))

    lt = data.get("locus_tag", {})
    if "prefix" not in lt:
        raise GffParseError(Diagnostic(Severity.ERROR, None, "config-missing",
                                       "config requires [locus_tag].prefix"))
    cds = data.get("cds", {})
    product = data.get("product", {})
    transcript = data.get("transcript", {})
    mode = transcript.get("mode", "nonredundant")
    if mode not in ("minimal", "nonredundant", "full"):
        raise GffParseError(Diagnostic(Severity.ERROR, None, "invalid-mode",
                                       f"transcript mode {mode!r} must be one of minimal/nonredundant/full"))

    cfg = MssConfig(
        locus_tag_prefix=lt["prefix"],
        locus_tag_width=lt.get("width", 6),
        locus_tag_start=lt.get("start", 10),
        locus_tag_step=lt.get("step", 10),
        transl_table=cds.get("transl_table", 1),
        product_default=product.get("default", "hypothetical protein"),
        product_map_path=product.get("map"),
        transcript_mode=mode,
        emit_mrna=transcript.get("emit_mrna", True),
    )
    return cfg, diags
