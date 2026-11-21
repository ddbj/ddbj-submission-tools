import pytest
from ddbj_submission_tools.nuc.ddbj_gff.errors import Severity, GffParseError
from ddbj_submission_tools.nuc.gff2mss.config import load_config


def _cfg(tmp_path, text):
    p = tmp_path / "c.toml"
    p.write_text(text)
    return load_config(str(p))


def test_load_config_defaults_and_unknown_key(tmp_path):
    cfg, diags = _cfg(tmp_path, '[locus_tag]\nprefix = "PFX"\n[unknown_section]\nx = 1\n')
    assert cfg.locus_tag_prefix == "PFX"
    assert cfg.locus_tag_width == 6 and cfg.locus_tag_start == 10 and cfg.locus_tag_step == 10
    assert cfg.transl_table == 1
    assert cfg.product_default == "hypothetical protein"
    assert [d.code for d in diags] == ["unknown-config-key"]
    assert diags[0].severity == Severity.WARNING


def test_missing_locus_tag_prefix_raises(tmp_path):
    with pytest.raises(GffParseError) as ei:
        _cfg(tmp_path, "[cds]\ntransl_table = 11\n")
    assert ei.value.diagnostic.code == "config-missing"


def test_metadata_sections_warn_moved_to_common(tmp_path):
    cfg, diags = _cfg(tmp_path, '[source]\norganism = "X"\n[assembly_gap]\nmin_length = 5\n'
                                '[locus_tag]\nprefix = "P"\n')
    assert [d.code for d in diags] == ["moved-to-common", "moved-to-common"]
    assert "SOURCE" in diags[0].message and "ASSEMBLY_GAP" in diags[1].message


def test_transcript_mode_default_and_explicit(tmp_path):
    cfg, _ = _cfg(tmp_path, '[locus_tag]\nprefix="P"\n')
    assert cfg.transcript_mode == "nonredundant"
    cfg, _ = _cfg(tmp_path, '[locus_tag]\nprefix="P"\n[transcript]\nmode="full"\n')
    assert cfg.transcript_mode == "full"


def test_invalid_transcript_mode_raises(tmp_path):
    with pytest.raises(GffParseError) as ei:
        _cfg(tmp_path, '[locus_tag]\nprefix="P"\n[transcript]\nmode="bogus"\n')
    assert ei.value.diagnostic.code == "invalid-mode"
