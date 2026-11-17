"""
Write DDBJ MSS annotation (.ann) files from a FASTA file only (no feature table).

Used by mss_builder. For assemblies that only need source + assembly_gap features.
"""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING, Optional

from common.common_builder import create_common
from common.fasta import parse_fasta_sequences
from common.gap_annotator import GapAnnotator, annotate_gaps
from common.submission_category import inject_defaults, validate_and_fill
from common.source_builder import (
    SequenceRoleEntry,
    ff_definition,
    source_qualifier,
    source_feature_rows,
)

if TYPE_CHECKING:
    from common.models import CommonModel

Row = list[str]


# ── Placeholder COMMON sections ───────────────────────────────────────────────

def _common_placeholder_wgs() -> list[Row]:
    """Placeholder COMMON rows for WGS mode (source included in COMMON)."""
    return [
        ["COMMON", "DBLINK",    "", "project",           ""],
        ["",       "",          "", "biosample",          ""],
        ["",       "SUBMITTER", "", "ab_name",            ""],
        ["",       "",          "", "contact",            ""],
        ["",       "",          "", "email",              ""],
        ["",       "",          "", "institute",          ""],
        ["",       "",          "", "country",            ""],
        ["",       "REFERENCE", "", "title",              ""],
        ["",       "",          "", "ab_name",            ""],
        ["",       "",          "", "status",             "Unpublished"],
        ["",       "",          "", "year",               ""],
        ["",       "DATE",      "", "hold_date",          ""],
        ["",       "source",    "1..E", "mol_type",       "genomic DNA"],
        ["",       "",          "", "ff_definition",      "@@[organism]@@ DNA, @@[submitter_seqid]@@"],
        ["",       "",          "", "submitter_seqid",    "@@[entry]@@"],
        ["",       "",          "", "organism",           ""],
    ]


def _common_placeholder_nonwgs() -> list[Row]:
    """Placeholder COMMON rows for non-WGS mode (source written per entry)."""
    return [
        ["COMMON", "DBLINK",    "", "project",     ""],
        ["",       "",          "", "biosample",    ""],
        ["",       "SUBMITTER", "", "ab_name",      ""],
        ["",       "",          "", "contact",      ""],
        ["",       "",          "", "email",        ""],
        ["",       "",          "", "institute",    ""],
        ["",       "",          "", "country",      ""],
        ["",       "REFERENCE", "", "title",        ""],
        ["",       "",          "", "ab_name",      ""],
        ["",       "",          "", "status",       "Unpublished"],
        ["",       "",          "", "year",         ""],
        ["",       "DATE",      "", "hold_date",    ""],
    ]


# ── Main writer ───────────────────────────────────────────────────────────────

def build_ann_rows(
    sequences: dict[str, str],
    common: Optional["CommonModel"] = None,
    sequence_roles: Optional[dict[str, SequenceRoleEntry]] = None,
    submission_category: Optional[str] = None,
    entry_rows: Optional[dict[str, list[Row]]] = None,
) -> list[Row]:
    """Rows of a DDBJ MSS .ann for *sequences* (id -> sequence, FASTA order).

    *entry_rows* maps an entry id to its feature rows (first column empty);
    they are placed after the entry's source (non-WGS) and before its gaps.
    See write_mss_ann for the WGS / non-WGS layout.
    """
    lengths = {seq_id: len(seq) for seq_id, seq in sequences.items()}
    all_ids = list(lengths.keys())

    gap_cfg = common.ASSEMBLY_GAP if common is not None else None
    gap_annotators: list[GapAnnotator] = []
    if gap_cfg:
        cfgs = gap_cfg if isinstance(gap_cfg, list) else [gap_cfg]
        gap_annotators = [
            GapAnnotator(
                linkage_evidence=cfg.linkage_evidence,
                min_gap_length=cfg.min_gap_length,
                max_gap_length=cfg.max_gap_length,
                gap_type=cfg.gap_type,
                estimated_length=cfg.estimated_length,
            )
            for cfg in cfgs
            if cfg.enabled
        ]
    gap_sequences = sequences if gap_annotators else {}

    # Determine WGS mode: no sequence_roles file, or every entry is unplaced
    def _is_unplaced(eid: str) -> bool:
        if sequence_roles is None:
            return True
        e = sequence_roles.get(eid)
        return e is None or e.type == "unplaced"

    is_wgs = all(_is_unplaced(eid) for eid in all_ids)

    chromosome_count = 0
    if sequence_roles:
        chromosome_count = sum(
            1 for e in sequence_roles.values() if e.type == "chromosome"
        )

    segment_count = 0
    if sequence_roles:
        segment_count = sum(
            1 for e in sequence_roles.values() if e.type == "segment"
        )

    # Base source qualifiers from common.SOURCE
    base_source: dict[str, str] = {}
    if common is not None and common.SOURCE:
        base_source.update(common.SOURCE)

    # Effective category: CLI --submission_category overrides JSON's _submission_category
    _extra = (common.model_extra or {}) if common is not None else {}
    effective_category: str = submission_category or _extra.get("_submission_category", "")

    source_id_key = common.SOURCE_IDENTIFIER if common is not None else None
    if not source_id_key and effective_category:
        from common.submission_category import get_category_rules
        source_id_key = get_category_rules(effective_category).source_identifier or None

    rows: list[Row] = []

    # ── COMMON section ────────────────────────────────────────────────────────
    if common is None:
        if is_wgs:
            rows.extend(_common_placeholder_wgs())
        else:
            rows.extend(_common_placeholder_nonwgs())
    else:
        common_dict = common.model_dump(exclude_none=True)
        if is_wgs:
            category = effective_category or "WGS"
            common_dict["_submission_category"] = category
            inject_defaults(common_dict, category)
            validate_and_fill(common_dict, category)
            rows.extend(create_common(common_dict, include_source=True))
        else:
            rows.extend(create_common(common_dict))

    # ── Per-entry body ────────────────────────────────────────────────────────
    for entry_id in all_ids:
        length = lengths[entry_id]
        location = f"1..{length}"

        role_entry: Optional[SequenceRoleEntry] = (
            sequence_roles.get(entry_id) if sequence_roles else None
        )
        is_circular = role_entry.is_circular if role_entry is not None else False

        if is_circular:
            rows.append([entry_id, "TOPOLOGY", "", "circular", ""])

        if not is_wgs:
            # Source feature per entry
            source_quals: dict[str, object] = dict(base_source)
            source_quals.update(source_qualifier(role_entry, entry_id, is_wgs=False, segment_count=segment_count))
            source_quals["ff_definition"] = ff_definition(
                role_entry, source_id_key, base_source.get("mol_type", ""),
                is_wgs=False, chromosome_count=chromosome_count, segment_count=segment_count,
            )

            source_entry_col = "" if is_circular else entry_id
            rows.extend(source_feature_rows(source_entry_col, location, source_quals))

        # Features (gff2mss), then assembly gaps. In WGS mode the source is in
        # COMMON, so the entry's first body row carries the entry name; in
        # non-WGS mode (or after a TOPOLOGY row) the first column stays empty.
        body: list[Row] = [list(r) for r in (entry_rows or {}).get(entry_id, [])]
        if gap_annotators and entry_id in gap_sequences:
            body.extend(annotate_gaps(gap_annotators, gap_sequences[entry_id], seq_name=None))
        if body and is_wgs and not is_circular:
            body[0][0] = entry_id
        rows.extend(body)

    return rows


def write_mss_ann(
    fsa_path: str,
    ann_path: str,
    common: Optional["CommonModel"] = None,
    sequence_roles: Optional[dict[str, SequenceRoleEntry]] = None,
    submission_category: Optional[str] = None,
) -> None:
    """
    Parse *fsa_path* (FASTA) and write a DDBJ MSS annotation file to *ann_path*.

    If *common* is provided (a validated CommonModel), its values are written
    into the COMMON section; otherwise placeholder lines are written.

    For WGS submissions (no *sequence_roles* provided, or all entries unplaced):
    - The source feature is placed in the COMMON block using ``@@[entry]@@`` and
      ``@@[organism]@@ DNA, @@[submitter_seqid]@@`` meta-notation.
    - Per-entry body contains only assembly_gap features (if any).

    For non-WGS submissions (*sequence_roles* provided with placed sequences):
    - Source features are written per entry with chromosome/organelle qualifiers.
    - Assembly_gap features follow the source for each entry.
    """
    sequences = parse_fasta_sequences(fsa_path)
    rows = build_ann_rows(sequences, common, sequence_roles, submission_category)

    with open(ann_path, "w") as fout:
        for row in rows:
            fout.write("\t".join(row) + "\n")

    print(f"[mss_builder] → {ann_path}", file=sys.stderr)
