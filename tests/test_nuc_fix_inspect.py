"""nuc.fix and nuc.inspect scripts on a small synthetic submission."""
import json
from pathlib import Path

from ddbj_submission_tools.nuc.fix import edit, empty_qualifier
from ddbj_submission_tools.nuc.inspect import cds_locus_tag, cds_seq, region, source_organism, submitter

FIX = Path(__file__).parent / "nuc_scripts_fixtures"
ANN, FASTA = FIX / "small.ann", FIX / "small.fasta"


def seq_of(fasta_text):
    return "".join(l for l in fasta_text.splitlines() if l and not l.startswith(">"))


def test_fix_edit_applies_json_rule(tmp_path, capsys):
    rules = tmp_path / "rules.json"
    rules.write_text(json.dumps({"rules": [{"action": "replace_value", "feature": "CDS", "qualifier": "product",
                                            "pattern": "example protein", "replacement": "example enzyme"}]}))
    assert edit.main([str(rules), str(ANN), "--outdir", str(tmp_path / "out")]) == 0
    text = (tmp_path / "out" / "small.ann").read_text()
    assert "\tproduct\texample enzyme" in text
    assert "example protein" not in text


def test_fix_empty_qualifier_dry_run_changes_nothing(tmp_path, capsys):
    ann = tmp_path / "e.ann"
    ann.write_text(ANN.read_text() + "\t5'UTR\t64..80\t\t\n")
    before = ann.read_text()
    assert empty_qualifier.main([str(ann), "--set", "*:note=tmp", "--dry-run"]) == 0
    assert ann.read_text() == before
    out = capsys.readouterr()
    assert "5'UTR" in out.out + out.err


def test_inspect_cds_seq_translates(tmp_path, capsys):
    out = tmp_path / "p.fa"
    cds_seq.main([str(ANN), "seq1", "--fasta", str(FASTA), "--translate", "-o", str(out)])
    first = out.read_text().split(">")[1].splitlines()[1]
    assert first.startswith("M")


def test_inspect_region_extracts_range(tmp_path, capsys):
    out = tmp_path / "r.fa"
    region.main([str(ANN), "seq1", "1..10", "--fasta", str(FASTA), "-o", str(out)])
    assert seq_of(out.read_text()).upper() == seq_of(FASTA.read_text().replace("//", ""))[:10].upper()


def test_inspect_cds_locus_tag_lists_cds(capsys):
    cds_locus_tag.main([str(ANN), "seq1", "--fasta", str(FASTA)])
    out = capsys.readouterr().out
    assert "1..63" in out
    assert "complement(join(100..150,200..250))" in out


def test_inspect_source_organism(capsys):
    assert source_organism.main([str(ANN), "seq1"]) == 0
    assert capsys.readouterr().out.strip() == "Examplea testa"


def test_inspect_submitter_email(capsys):
    assert submitter.main([str(ANN), "--field", "email"]) == 0
    assert capsys.readouterr().out.strip() == "example@example.org"


def test_inspect_submitter_ab_name(capsys):
    assert submitter.main([str(ANN), "--field", "ab_name"]) == 0
    assert capsys.readouterr().out.strip() == "Mishima,T."
