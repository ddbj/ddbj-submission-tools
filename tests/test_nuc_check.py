"""nuc.check scripts on a small synthetic submission."""
import subprocess
import sys
from pathlib import Path

from ddbj_submission_tools.nuc.check import ann_set, mat_peptide, review, rrna, short_complete_cds, trna

FIX = Path(__file__).parent / "nuc_scripts_fixtures"
ANN, FASTA = FIX / "small.ann", FIX / "small.fasta"


def variant(tmp_path, old, new, name="v.ann"):
    p = tmp_path / name
    text = ANN.read_text()
    assert old in text
    p.write_text(text.replace(old, new))
    return p


def test_mat_peptide_in_cds_is_fine(capsys):
    assert mat_peptide.main([str(ANN), "seq1", "--fasta", str(FASTA)]) == 0


def test_mat_peptide_outside_cds_reported(tmp_path, capsys):
    ann = variant(tmp_path, "mat_peptide\t4..30", "mat_peptide\t70..96")
    assert mat_peptide.main([str(ann), "seq1", "--fasta", str(FASTA)]) == 1
    assert "Error: line 21:" in capsys.readouterr().err


def test_short_complete_cds_lists_20aa(capsys):
    assert short_complete_cds.main([str(ANN), "--min-aa", "30"]) == 0
    rows = [l for l in capsys.readouterr().out.splitlines() if l and not l.startswith("#")]
    assert len(rows) == 1                     # the 20 aa CDS only (header line starts with '#')
    assert rows[0].split("\t")[:4] == ["seq1", "1..63", "63", "20"]


def test_review_reports_spaced_qualifier(tmp_path, capsys):
    ann = variant(tmp_path, "\t\t\tstrain\tJCM1234", "\t\t\tstrain name\tJCM1234")
    review.main([str(ann)])
    out = capsys.readouterr().out
    assert "qualifier name with space" in out
    assert "strain name" in out


def test_ann_set_reports_differing_common(tmp_path, capsys):
    (tmp_path / "a.ann").write_text(ANN.read_text())
    variant(tmp_path, "PRJDB00001", "PRJDB00002", "b.ann")
    assert ann_set.main([str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert any(l.split()[:3] == ["DIFFER", "DBLINK", "project"] for l in out.splitlines())


def test_rrna_extracts_fasta(tmp_path, capsys):
    out = tmp_path / "r.fa"
    rrna.main([str(ANN), "seq1", "--fasta", str(FASTA), "-o", str(out)])
    assert out.read_text().count(">") == 1


def test_trna_extracts_fasta(tmp_path, capsys):
    out = tmp_path / "t.fa"
    trna.main([str(ANN), "seq1", "--fasta", str(FASTA), "-o", str(out)])
    assert out.read_text().count(">") == 1


def test_python_m_runs():
    r = subprocess.run([sys.executable, "-m", "ddbj_submission_tools.nuc.check.review", str(ANN)],
                       capture_output=True, text=True)
    assert r.returncode in (0, 1), r.stderr
    assert "== Findings" in r.stdout
