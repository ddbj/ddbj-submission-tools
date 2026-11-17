"""The gff2mss example in examples/gff2mss runs as the README shows."""
from pathlib import Path

from gff2mss.cli import main

EX = Path(__file__).resolve().parents[1] / "examples" / "gff2mss"


def test_gff2mss_example_runs(tmp_path):
    rc = main(["--gff", str(EX / "citrus.gff3"), "--fasta", str(EX / "citrus.fasta"),
               "--mss-config", str(EX / "mss.toml"), "--common", str(EX / "common.json"),
               "--out", str(tmp_path / "citrus")])
    assert rc == 0
    ann = (tmp_path / "citrus.ann").read_text()
    assert ann.startswith("COMMON")
    assert "\tCDS\t" in ann
