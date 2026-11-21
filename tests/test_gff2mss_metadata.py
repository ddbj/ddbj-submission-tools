"""gff2mss writes metadata exactly like mss_builder; its config holds annotation settings only."""
import sys
from pathlib import Path

from ddbj_submission_tools.nuc.gff2mss.cli import main as gff2mss_main
from ddbj_submission_tools.nuc.mss_builder.cli import main as mss_builder_main

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


def test_mss_style_fasta_reads_like_mss_builder(tmp_path, monkeypatch):
    # MSS FASTA (records end with "//") and lcl| headers: same parsing as mss_builder
    src = (EX / "example.gnm.fasta").read_text().replace(">", ">lcl|")
    mss_fasta = src.replace("\n>", "\n//\n>").rstrip("\n") + "\n//\n"
    fasta = _write(tmp_path, "mss.fasta", mss_fasta)
    argv = ["--gff", _write(tmp_path, "in.gff3", EMPTY_GFF), "--fasta", fasta,
            "--mss-config", _write(tmp_path, "mss.toml", MIN_TOML),
            "--common", str(EX / "example.gnm.common.json"), "--sequence_roles", ROLES,
            "--out", str(tmp_path / "g2m")]
    assert gff2mss_main(argv) == 0
    monkeypatch.setattr(sys, "argv", ["mss_builder", "--common", str(EX / "example.gnm.common.json"),
                                      "--sequence_roles", ROLES, "-o", str(tmp_path / "mb"), "-p", "gnm", fasta])
    mss_builder_main()
    assert (tmp_path / "g2m.ann").read_text() == (tmp_path / "mb" / "gnm.ann").read_text()
    assert "/" not in "".join(l for l in (tmp_path / "g2m.fasta").read_text().splitlines()
                              if not l.startswith(">") and l != "//")


def _gff2mss_exit(tmp_path, extra=(), common=None):
    argv = ["--gff", _write(tmp_path, "in.gff3", EMPTY_GFF), "--fasta", str(EX / "example.gnm.fasta"),
            "--mss-config", _write(tmp_path, "mss.toml", MIN_TOML),
            "--common", common or str(EX / "example.gnm.common.json"), "--out", str(tmp_path / "g2m"), *extra]
    try:
        gff2mss_main(argv)
    except SystemExit as exc:
        return exc.code
    raise AssertionError("gff2mss did not exit")


def test_malformed_common_json_is_a_clean_error(tmp_path):
    msg = _gff2mss_exit(tmp_path, common=_write(tmp_path, "bad.json", '{"DBLINK": '))
    assert isinstance(msg, str) and "failed to load --common file" in msg and "bad.json" in msg


def test_invalid_common_field_is_a_clean_error(tmp_path):
    msg = _gff2mss_exit(tmp_path, common=_write(tmp_path, "bad.json", '{"SUBMITTER": "not an object"}'))
    assert isinstance(msg, str) and "failed to load --common file" in msg


def test_malformed_sequence_roles_is_a_clean_error(tmp_path):
    roles = _write(tmp_path, "roles.tsv", "only_one_column\n")
    msg = _gff2mss_exit(tmp_path, extra=["--sequence_roles", roles])
    assert isinstance(msg, str) and "failed to load sequence role file" in msg and "roles.tsv" in msg


def test_unannotated_wgs_sequences_are_in_output_fasta(tmp_path):
    from ddbj_submission_tools.nuc.common.fasta import parse_fasta_sequences
    assert _gff2mss(tmp_path, "wgs") == 0
    out = parse_fasta_sequences(str(tmp_path / "g2m.fasta"))
    inp = parse_fasta_sequences(str(EX / "example.wgs.fasta"))
    assert list(out) == list(inp)
    assert {k: v.upper() for k, v in out.items()} == {k: v.upper() for k, v in inp.items()}
