from pathlib import Path

from ddbj_submission_tools.nuc.ddbj_gff import parse
from ddbj_submission_tools.nuc.ddbj_gff.tools._base import run_tool

FX = Path(__file__).parent / "fixtures" / "tools"


def _run(mod, fixture, tmp_path):
    out, rep = tmp_path / "o.gff3", tmp_path / "r.txt"
    assert run_tool(mod.SPEC, [str(FX / fixture), "--out", str(out), "--report", str(rep)]) == 0
    return parse(out.read_text()), rep.read_text()


def test_liftoff_cds_get_phases(tmp_path):
    from ddbj_submission_tools.nuc.ddbj_gff.tools import liftoff
    doc, rep = _run(liftoff, "liftoff.gff3", tmp_path)
    cds = [f for f in doc.features if f.type == "CDS"]
    assert cds and all(s.phase in (0, 1, 2) for f in cds for s in f.spans)
    assert rep.count("set-phase") == len(cds)
    assert [r.id for r in doc.roots] == ["Mp1g00070"]


def test_maker_drops_evidence_features(tmp_path):
    from ddbj_submission_tools.nuc.ddbj_gff.tools import maker
    doc, rep = _run(maker, "maker.gff3", tmp_path)
    assert {f.type for f in doc.features} == {"gene", "mRNA", "exon", "CDS"}
    assert all(r.type == "gene" for r in doc.roots)
    for t in ("match", "match_part", "contig"):
        assert f"drop-feature {t}" in rep


def test_snap_hierarchy_with_phases(tmp_path):
    from ddbj_submission_tools.nuc.ddbj_gff.tools import snap
    doc, rep = _run(snap, "snap.gff3", tmp_path)
    assert all(r.type == "gene" for r in doc.roots) and len(doc.roots) == 9
    cds = doc.feature_index["chr1-snap.3.t1.cds"]
    assert len(cds.spans) == 2 and all(s.phase is not None for s in cds.spans)
    assert "set-phase" in rep
