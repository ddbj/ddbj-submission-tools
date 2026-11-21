from pathlib import Path

from ddbj_submission_tools.nuc.ddbj_gff.tools.cli import TOOLS, main

FIX = Path(__file__).parent / "fixtures" / "tools"


def test_lists_only_tool_modules():
    assert "augustus" in TOOLS and "transdecoder" in TOOLS
    assert not {"prepass", "snap_reader", "cli"} & set(TOOLS)


def test_dispatches_to_tool(tmp_path):
    out = tmp_path / "a.gff3"
    assert main(["augustus", str(FIX / "augustus_stopExcluded.gff3"), "--out", str(out),
                 "--report", str(tmp_path / "r.txt")]) == 0
    assert "\tmRNA\t" in out.read_text()


def test_unknown_tool_is_usage_error(capsys):
    assert main(["nosuch", "x"]) == 2
    assert "tools:" in capsys.readouterr().err
