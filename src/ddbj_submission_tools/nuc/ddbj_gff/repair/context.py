from __future__ import annotations

from dataclasses import dataclass


@dataclass
class RepairContext:
    sequences: dict | None = None   # seqid -> Bio.Seq.Seq (nucleotide)
    transl_table: int = 1           # default table when a CDS omits transl_table
    inference_value: str = "ab initio prediction"   # /inference value written verbatim
    inference_types: str = "mRNA"   # comma-separated feature types; a "*" entry matches any
    inference_sources: str = ""     # comma-separated GFF column-2 sources; empty = any
    min_cds_length: int = 90        # a shorter total CDS is a short-CDS candidate (ANN4300)
