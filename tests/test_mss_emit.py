from ddbj_submission_tools.nuc.gff2mss.model import MssQualifier, MssFeature
from ddbj_submission_tools.nuc.gff2mss.emit import emit_fasta, feature_rows


def test_emit_fasta_wrap_and_terminator():
    seqs = {"chr1": "A" * 130}
    out = emit_fasta(seqs)
    lines = out.splitlines()
    assert lines[0] == ">chr1"
    assert lines[1] == "A" * 60
    assert lines[2] == "A" * 60
    assert lines[3] == "A" * 10
    assert lines[4] == "//"


def test_emit_fasta_multi_entry():
    out = emit_fasta({"chr1": "AAA", "chr2": "CCC"})
    assert out.splitlines() == [">chr1", "AAA", "//", ">chr2", "CCC", "//"]


def test_feature_rows_first_row_carries_key_and_location():
    feat = MssFeature("CDS", "1..9", [MssQualifier("locus_tag", "L_1"),
                                       MssQualifier("product", "x")])
    rows = feature_rows(feat)
    assert rows[0] == ["", "CDS", "1..9", "locus_tag", "L_1"]
    assert rows[1] == ["", "", "", "product", "x"]
