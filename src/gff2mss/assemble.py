"""Assemble a DDBJ MSS .ann: features from gff2mss, metadata exactly as mss_builder writes it."""
from __future__ import annotations

import sys

from Bio.Seq import Seq

from ddbj_gff import parse
from ddbj_gff.io import open_text
from gff2mss.config import load_config
from gff2mss.convert import build_entry_features
from gff2mss.emit import feature_rows
from gff2mss.product_map import load_product_map

from common.fasta import parse_fasta_sequences
from common.models import load_common_json
from common.source_builder import load_sequence_roles
from mss_builder.ann_writer import build_ann_rows


def build_ann_text(gff_path, fasta_path, mss_config_path, common_path,
                   sequence_roles_path, submission_category, locus_tag_start=None):
    cfg, config_diags = load_config(mss_config_path)
    for d in config_diags:
        sys.stderr.write(f"[gff2mss] {d.severity.value} {d.code}: {d.message}\n")
    if cfg.product_map_path:
        cfg.product_map = load_product_map(cfg.product_map_path)
    if locus_tag_start is not None:
        cfg.locus_tag_start = locus_tag_start  # continue numbering across compartments (genome-unique)

    # same FASTA reading as mss_builder (MSS "//" separators, lcl| prefixes)
    seq_text = parse_fasta_sequences(fasta_path)
    seqs = {k: Seq(v) for k, v in seq_text.items()}
    with open_text(gff_path) as fh:
        doc = parse(fh.read())

    diagnostics: list = []
    per_entry = build_entry_features(doc, seqs, cfg, diagnostics)
    for d in diagnostics:
        if d.code in ("root-ignored", "missing-sequence"):      # features silently lost otherwise
            sys.stderr.write(f"[gff2mss] {d.severity.value} {d.code}: {d.message}\n")

    entry_rows = {eid: [row for feat in feats for row in feature_rows(feat)]
                  for eid, feats in per_entry.items()}
    rows = build_ann_rows(seq_text, load_common_json(common_path),
                          load_sequence_roles(sequence_roles_path) if sequence_roles_path else None,
                          submission_category or None, entry_rows)
    ann_text = "\n".join("\t".join(r) for r in rows) + "\n"
    # output every input sequence (annotation or not), in FASTA order
    return ann_text, seq_text
