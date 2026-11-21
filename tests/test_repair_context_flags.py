import os
from ddbj_submission_tools.nuc.ddbj_gff.repair.context import RepairContext
from ddbj_submission_tools.nuc.ddbj_gff.repair.cli import main

FIX = os.path.join(os.path.dirname(__file__), "fixtures")
GFF = os.path.join(FIX, "repair_internal_stop.gff3")


def test_context_defaults():
    ctx = RepairContext()
    assert ctx.inference_value == "ab initio prediction"
    assert ctx.inference_types == "mRNA"
    assert ctx.inference_sources == ""
    assert ctx.min_cds_length == 90


def test_context_keeps_its_existing_fields_and_positional_order():
    # existing call sites construct RepairContext() and RepairContext(sequences=..., ...)
    ctx = RepairContext(sequences={"c1": "AAA"}, transl_table=11)
    assert ctx.sequences == {"c1": "AAA"} and ctx.transl_table == 11
    assert ctx.min_cds_length == 90          # new fields still defaulted


def test_cli_accepts_the_four_new_flags(tmp_path, capsys):
    out = tmp_path / "out.gff3"
    rc = main(["--gff", GFF, "--apply", "utr-absent-to-partial-mrna",
               "--inference-value", "similar to RNA sequence:INSD:LC123456.1",
               "--inference-types", "mRNA,tRNA",
               "--inference-sources", "AUGUSTUS",
               "--min-cds-length", "30",
               "--out", str(out), "--report", str(tmp_path / "r.txt")])
    assert rc == 0
    assert out.exists()


def test_cli_passes_the_flags_into_the_context(tmp_path, monkeypatch):
    seen = {}
    import ddbj_submission_tools.nuc.ddbj_gff.repair.cli as cli_mod

    real = cli_mod.RepairContext

    def spy(**kw):
        seen.update(kw)
        return real(**kw)

    monkeypatch.setattr(cli_mod, "RepairContext", spy)
    main(["--gff", GFF, "--apply", "utr-absent-to-partial-mrna",
          "--inference-value", "EXISTENCE:none", "--inference-types", "*",
          "--inference-sources", "BRAKER,AUGUSTUS", "--min-cds-length", "42",
          "--out", str(tmp_path / "o.gff3"), "--report", str(tmp_path / "r.txt")])
    assert seen["inference_value"] == "EXISTENCE:none"
    assert seen["inference_types"] == "*"
    assert seen["inference_sources"] == "BRAKER,AUGUSTUS"
    assert seen["min_cds_length"] == 42


def test_cli_defaults_match_the_context_defaults(tmp_path, monkeypatch):
    seen = {}
    import ddbj_submission_tools.nuc.ddbj_gff.repair.cli as cli_mod
    real = cli_mod.RepairContext

    def spy(**kw):
        seen.update(kw)
        return real(**kw)

    monkeypatch.setattr(cli_mod, "RepairContext", spy)
    main(["--gff", GFF, "--apply", "utr-absent-to-partial-mrna",
          "--out", str(tmp_path / "o.gff3"), "--report", str(tmp_path / "r.txt")])
    d = RepairContext()
    assert seen["inference_value"] == d.inference_value
    assert seen["inference_types"] == d.inference_types
    assert seen["inference_sources"] == d.inference_sources
    assert seen["min_cds_length"] == d.min_cds_length
