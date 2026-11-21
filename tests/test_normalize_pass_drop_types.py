from ddbj_submission_tools.nuc.ddbj_gff import parse
from ddbj_submission_tools.nuc.ddbj_gff.normalize.config import NormalizeConfig, load_normalize_config
from ddbj_submission_tools.nuc.ddbj_gff.normalize.normalize import normalize

GFF = """##gff-version 3
chr1\tA\tgene\t1\t20\t.\t+\t.\tID=g1
chr1\tA\tmRNA\t1\t20\t.\t+\t.\tID=g1.t1;Parent=g1
chr1\tA\tCDS\t1\t9\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1
chr1\tA\tCDS\t15\t20\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1
chr1\tA\tintron\t10\t14\t.\t+\t.\tID=g1.t1.i1;Parent=g1.t1
chr1\tA\tstop_codon\t18\t20\t.\t+\t0\tID=g1.t1.stop;Parent=g1.t1
chr1\tA\tmatch\t30\t40\t.\t+\t.\tID=m1
chr1\tA\tmatch_part\t30\t40\t.\t+\t.\tID=m1.p1;Parent=m1
"""


def test_drops_types_and_reports():
    norm, rep = normalize(parse(GFF), config=NormalizeConfig(taxid=1, drop_types=("intron", "stop_codon", "match")))
    assert {f.type for f in norm.features} == {"gene", "mRNA", "CDS"}
    assert "g1.t1.i1" not in norm.feature_index and "m1.p1" not in norm.feature_index
    assert [c.type for c in norm.feature_index["g1.t1"].children] == ["CDS"]
    assert [r.id for r in norm.roots] == ["g1"]
    msgs = {(c.action, c.target) for c in rep.applied}
    assert ("drop-feature", "intron") in msgs and ("drop-feature", "stop_codon") in msgs
    assert ("drop-feature", "match") in msgs


def test_default_off_keeps_everything():
    norm, _ = normalize(parse(GFF), config=NormalizeConfig(taxid=1))
    assert "g1.t1.i1" in norm.feature_index


def test_toml_key(tmp_path):
    p = tmp_path / "n.toml"
    p.write_text('[normalize]\ndrop_types = ["intron"]\n')
    assert load_normalize_config(str(p)).drop_types == ("intron",)
