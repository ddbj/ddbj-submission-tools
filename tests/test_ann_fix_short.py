from ann_fix.annfile import AnnFile
from ann_fix.cli import main
from ann_fix.short import drop_short, short_report


def _rows(*rows):
    return "".join("\t".join(r) + "\n" for r in rows)


def _files(tmp_path, with_entries=True):
    ann = tmp_path / "x.ann"
    common = [("COMMON", "DATE", "", "hold_date", "20270101")]
    entries = [("seq1", "source", "1..1200", "organism", "E. coli"), ("", "CDS", "1..900", "product", "p1"),
               ("", "tRNA", "1000..1075", "product", "tRNA-Ala"),
               ("seq2", "source", "1..500", "organism", "E. coli"), ("", "tRNA", "10..85", "product", "tRNA-Ile"),
               ("", "CDS", "100..399", "product", "p2"),
               ("seq3", "source", "1..800", "organism", "E. coli"), ("", "tRNA", "10..85", "product", "tRNA-Ala")]
    ann.write_text(_rows(*common, *(entries if with_entries else [])))
    fa = tmp_path / "x.fasta"
    fa.write_text(">seq1\n" + "A" * 60 + "\n" + "A" * 1140 + "\n//\n>seq2\n" + "C" * 500 + "\n//\n>seq3\n" + "G" * 800 + "\n//\n")
    return ann, fa


def test_short_report_counts_what_would_be_lost(tmp_path):
    ann, fa = _files(tmp_path)
    r = short_report(AnnFile.read(ann), fa, 1000)
    assert r["short"] == [{"entry": "seq2", "length": 500}, {"entry": "seq3", "length": 800}]
    assert r["count"] == 2 and r["sequences"] == 3 and r["short_bases"] == 1300 and r["total_bases"] == 2500
    assert r["lost_features"] == {"tRNA": 2, "CDS": 1}
    assert r["lost_trna_products"] == ["tRNA-Ile"]          # tRNA-Ala still has a copy on seq1


def test_drop_short_removes_entries_from_ann_and_fasta(tmp_path):
    ann, fa = _files(tmp_path)
    a = AnnFile.read(ann)
    fasta_text = drop_short(a, fa, 1000)
    assert [f.entry for f in a.features()] == ["COMMON", "seq1", "seq1", "seq1"]
    assert fasta_text == ">seq1\n" + "A" * 60 + "\n" + "A" * 1140 + "\n//\n"


def test_drop_short_with_common_source_only(tmp_path):
    ann, fa = _files(tmp_path, with_entries=False)
    a = AnnFile.read(ann)
    assert drop_short(a, fa, 1000).count(">") == 1
    assert [f.key for f in a.features()] == ["DATE"]


def test_cli_short_report_and_drop(tmp_path, capsys):
    ann, fa = _files(tmp_path)
    assert main(["--ann", str(ann), "--fasta", str(fa), "--short-seqs"]) == 0
    out = capsys.readouterr().out
    assert "2 of 3 sequences" in out and "tRNA-Ile" in out
    o_ann, o_fa = tmp_path / "o.ann", tmp_path / "o.fasta"
    assert main(["--ann", str(ann), "--fasta", str(fa), "--drop-short", "--out", str(o_ann), "--fasta-out", str(o_fa)]) == 0
    assert o_fa.read_text().count(">") == 1 and "seq2" not in o_ann.read_text()
    assert main(["--ann", str(ann), "--fasta", str(fa), "--drop-short", "--out", str(o_ann)]) == 2   # needs --fasta-out
