"""gff2mss writes metadata exactly like mss_builder; its config holds annotation settings only."""
import sys
from pathlib import Path

from gff2mss.cli import main as gff2mss_main
from mss_builder.cli import main as mss_builder_main

EX = Path(__file__).resolve().parents[1] / "examples" / "mss_builder"
ROLES = str(EX / "example.sequence_roles.tsv")
EMPTY_GFF = "##gff-version 3\n"
MIN_TOML = '[locus_tag]\nprefix="TEST"\n'


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def _gff2mss(tmp_path, kind, toml=MIN_TOML, gff=EMPTY_GFF, extra=()):
    argv = ["--gff", _write(tmp_path, "in.gff3", gff), "--fasta", str(EX / f"example.{kind}.fasta"),
            "--mss-config", _write(tmp_path, "mss.toml", toml),
            "--common", str(EX / f"example.{kind}.common.json"), "--out", str(tmp_path / "g2m"), *extra]
    return gff2mss_main(argv)


def _same_as_mss_builder(tmp_path, monkeypatch, kind, g2m_extra, mb_extra):
    assert _gff2mss(tmp_path, kind, extra=g2m_extra) == 0
    monkeypatch.setattr(sys, "argv", ["mss_builder", "--common", str(EX / f"example.{kind}.common.json"),
                                      *mb_extra, "-o", str(tmp_path / "mb"), "-p", kind,
                                      str(EX / f"example.{kind}.fasta")])
    mss_builder_main()
    assert (tmp_path / "g2m.ann").read_text() == (tmp_path / "mb" / f"{kind}.ann").read_text()


def test_empty_gff_matches_mss_builder_wgs(tmp_path, monkeypatch):
    _same_as_mss_builder(tmp_path, monkeypatch, "wgs", [], [])


def test_empty_gff_matches_mss_builder_gnm(tmp_path, monkeypatch):
    _same_as_mss_builder(tmp_path, monkeypatch, "gnm", ["--sequence-roles", ROLES], ["--sequence_roles", ROLES])


def test_old_source_section_warns_and_converts(tmp_path, capsys):
    toml = MIN_TOML + '[source]\norganism="x"\n[assembly_gap]\nmin_length=10\n'
    assert _gff2mss(tmp_path, "wgs", toml=toml) == 0
    assert capsys.readouterr().err.count("moved-to-common") == 2


def test_cli_accepts_underscore_aliases(tmp_path):
    assert _gff2mss(tmp_path, "gnm", extra=["--sequence_roles", ROLES, "--submission_category", "GNM"]) == 0
    assert "chromosome" in (tmp_path / "g2m.ann").read_text()


def test_missing_sequence_still_reported(tmp_path, capsys):
    gff = EMPTY_GFF + "nosuchseq\tx\tgene\t1\t90\t.\t+\t.\tID=g1\n"
    assert _gff2mss(tmp_path, "wgs", gff=gff) == 0
    assert "missing-sequence" in capsys.readouterr().err


def test_out_directory_is_created(tmp_path):
    argv = ["--gff", _write(tmp_path, "in.gff3", EMPTY_GFF), "--fasta", str(EX / "example.wgs.fasta"),
            "--mss-config", _write(tmp_path, "mss.toml", MIN_TOML),
            "--common", str(EX / "example.wgs.common.json"), "--out", str(tmp_path / "new" / "dir" / "x")]
    assert gff2mss_main(argv) == 0
    assert (tmp_path / "new" / "dir" / "x.ann").exists()
