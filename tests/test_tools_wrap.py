import importlib
from pathlib import Path

import pytest

from ddbj_submission_tools.nuc.ddbj_gff import parse
from ddbj_submission_tools.nuc.ddbj_gff.tools._base import run_tool

FX = Path(__file__).parent / "fixtures" / "tools"


@pytest.mark.parametrize("tool,fixture", [("glimmerhmm", "glimmerhmm.gff"), ("ginger", "ginger.gff")])
def test_parentless_mrna_gets_gene(tool, fixture, tmp_path):
    mod = importlib.import_module(f"ddbj_submission_tools.nuc.ddbj_gff.tools.{tool}")
    out, rep = tmp_path / "o.gff3", tmp_path / "r.txt"
    assert run_tool(mod.SPEC, [str(FX / fixture), "--out", str(out), "--report", str(rep)]) == 0
    doc = parse(out.read_text())
    n_in = sum(1 for l in (FX / fixture).read_text().splitlines() if "\tmRNA\t" in l)
    assert n_in >= 2
    assert all(r.type == "gene" for r in doc.roots) and len(doc.roots) == n_in
    mrnas = [f for f in doc.features if f.type == "mRNA"]
    assert len(mrnas) == n_in and all(m.parents and m.parents[0].type == "gene" for m in mrnas)
    assert rep.read_text().count("add-parent-feature") == n_in
