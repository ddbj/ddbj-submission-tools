from pathlib import Path
from ddbj_submission_tools.nuc.gff2mss.cli import main

FIX = Path(__file__).parent / "mss_fixtures"


def test_snapshot_ann_and_fasta(tmp_path):
    out = tmp_path / "got"
    rc = main(["--gff", str(FIX / "mini.gff3"), "--fasta", str(FIX / "mini.fa"),
               "--mss-config", str(FIX / "config.toml"), "--common", str(FIX / "common.json"),
               "--sequence_roles", str(FIX / "sequence_roles.tsv"), "--out", str(out)])
    assert rc == 0
    assert (tmp_path / "got.ann").read_text() == (FIX / "expected.ann").read_text()
    assert (tmp_path / "got.fasta").read_text() == (FIX / "expected.fasta").read_text()
