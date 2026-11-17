"""build_ann_rows: the metadata rows mss_builder writes, with optional per-entry feature rows."""
from pathlib import Path

from common.fasta import parse_fasta_sequences
from common.models import load_common_json
from common.source_builder import load_sequence_roles
from mss_builder.ann_writer import build_ann_rows, write_mss_ann

EX = Path(__file__).resolve().parents[1] / "examples" / "mss_builder"
ROLES = load_sequence_roles(str(EX / "example.sequence_roles.tsv"))


def _inputs(kind):
    fasta = EX / f"example.{kind}.fasta"
    return fasta, parse_fasta_sequences(str(fasta)), load_common_json(str(EX / f"example.{kind}.common.json"))


def _text(rows):
    return "".join("\t".join(r) + "\n" for r in rows)


def _entry_block(rows, entry_id):
    """Rows from the first row naming entry_id up to the next row naming another entry."""
    start = next(i for i, r in enumerate(rows) if r[0] == entry_id)
    end = next((i for i in range(start + 1, len(rows)) if rows[i][0] not in ("", entry_id)), len(rows))
    return rows[start:end]


def test_no_features_matches_write_mss_ann_wgs(tmp_path):
    fasta, seqs, common = _inputs("wgs")
    write_mss_ann(str(fasta), str(tmp_path / "a.ann"), common=common)
    assert _text(build_ann_rows(seqs, common, None, None)) == (tmp_path / "a.ann").read_text()


def test_no_features_matches_write_mss_ann_gnm(tmp_path):
    fasta, seqs, common = _inputs("gnm")
    write_mss_ann(str(fasta), str(tmp_path / "a.ann"), common=common, sequence_roles=ROLES)
    assert _text(build_ann_rows(seqs, common, ROLES, None)) == (tmp_path / "a.ann").read_text()


def test_wgs_feature_row_carries_entry_name():
    _, seqs, common = _inputs("wgs")
    feat = ["", "gene", "1..90", "locus_tag", "T_000010"]
    rows = build_ann_rows(seqs, common, None, None, entry_rows={"contig01": [feat]})
    block = _entry_block(rows, "contig01")
    assert block[0] == ["contig01", "gene", "1..90", "locus_tag", "T_000010"]
    gaps = [r for r in block if r[1] == "assembly_gap"]
    assert gaps and all(r[0] == "" for r in gaps)
    assert sum(r[0] == "contig01" for r in rows) == 1


def test_wgs_entry_without_features_or_gaps_has_no_rows():
    _, seqs, common = _inputs("wgs")
    no_gap = {k: v for k, v in seqs.items() if "n" not in v.lower()}
    assert no_gap, "example needs a contig without gaps"
    rows = build_ann_rows(no_gap, common, None, None)
    assert not any(r[0] in no_gap for r in rows)


def test_circular_entry_with_features():
    _, seqs, common = _inputs("gnm")
    feat = ["", "gene", "1..90", "locus_tag", "T_000010"]
    rows = build_ann_rows(seqs, common, ROLES, None, entry_rows={"mt": [feat]})
    i = rows.index(["mt", "TOPOLOGY", "", "circular", ""])
    assert rows[i + 1][:2] == ["", "source"]
    src_end = next(j for j in range(i + 1, len(rows)) if rows[j][1] not in ("", "source"))
    assert rows[src_end] == feat
