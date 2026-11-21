"""Test-only: render gff2mss features as .ann text with a bare source per entry (no COMMON)."""
from ddbj_submission_tools.nuc.gff2mss.convert import build_entry_features
from ddbj_submission_tools.nuc.gff2mss.emit import feature_rows


def ann_text(doc, seqs, cfg):
    rows = []
    for seqid, feats in build_entry_features(doc, seqs, cfg, []).items():
        rows.append([seqid, "source", f"1..{len(seqs[seqid])}", "submitter_seqid", seqid])
        rows.extend(row for feat in feats for row in feature_rows(feat))
    return "".join("\t".join(r) + "\n" for r in rows)
