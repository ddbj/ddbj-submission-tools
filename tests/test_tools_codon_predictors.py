import importlib
from pathlib import Path

import pytest

from ddbj_submission_tools.nuc.ddbj_gff import parse
from ddbj_submission_tools.nuc.ddbj_gff.tools._base import run_tool

FX = Path(__file__).parent / "fixtures" / "tools"
CASES = [("augustus", "augustus_stopExcluded.gff3"), ("braker", "braker3.gtf"),
         ("genemark", "genemark_et.gtf"), ("tiberius", "tiberius.gtf")]


@pytest.mark.parametrize("tool,fixture", CASES)
def test_codon_predictor_outputs_gene_mrna_cds(tool, fixture, tmp_path):
    mod = importlib.import_module(f"ddbj_submission_tools.nuc.ddbj_gff.tools.{tool}")
    out, rep = tmp_path / "o.gff3", tmp_path / "r.txt"
    assert run_tool(mod.SPEC, [str(FX / fixture), "--out", str(out), "--report", str(rep)]) == 0
    doc = parse(out.read_text())
    assert {f.type for f in doc.features} <= {"gene", "mRNA", "CDS", "exon"}
    mrnas = [f for f in doc.features if f.type == "mRNA"]
    assert mrnas and all(sum(1 for c in m.children for _ in [c] if c.type == "CDS") == 1 for m in mrnas)
    assert all(r.type == "gene" for r in doc.roots)
    assert "drop-feature stop_codon" in rep.read_text()


def test_augustus_stop_excluded_cds_is_extended(tmp_path):
    from ddbj_submission_tools.nuc.ddbj_gff.tools import augustus
    out, rep = tmp_path / "o.gff3", tmp_path / "r.txt"
    run_tool(augustus.SPEC, [str(FX / "augustus_stopExcluded.gff3"), "--out", str(out), "--report", str(rep)])
    cds = parse(out.read_text()).feature_index["g1.t1.cds"]
    assert min(s.start for s in cds.spans) == 1807        # stop_codon 1807..1809 folded in (minus strand)
    assert "extend-cds g1.t1.cds" in rep.read_text()
