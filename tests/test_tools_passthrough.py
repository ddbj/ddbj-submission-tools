import importlib
from pathlib import Path

import pytest

from ddbj_gff import parse, write
from ddbj_gff.tools._base import run_tool

FX = Path(__file__).parent / "fixtures" / "tools"
CASES = [("annevo", "annevo.gff3"), ("helixer", "helixer.gff3"), ("evm", "evm.gff3"),
         ("funannotate", "funannotate.gff3"), ("transdecoder", "transdecoder.gff3"),
         ("pasa", "pasa.gff3"), ("gemoma", "gemoma.gff")]


@pytest.mark.parametrize("tool,fixture", CASES)
def test_passthrough_is_already_canonical(tool, fixture, tmp_path):
    mod = importlib.import_module(f"ddbj_gff.tools.{tool}")
    out = tmp_path / "o.gff3"
    assert run_tool(mod.SPEC, [str(FX / fixture), "--out", str(out)]) == 0
    doc = parse(out.read_text())
    assert doc.roots and all(r.type == "gene" for r in doc.roots)
    mrnas = [f for f in doc.features if f.type == "mRNA"]
    assert mrnas and all(any(c.type == "CDS" for c in m.children) for m in mrnas)
    src = parse((FX / fixture).read_text())
    assert out.read_text() == write(src) or out.read_text() == "##gff-version 3\n" + write(src)
