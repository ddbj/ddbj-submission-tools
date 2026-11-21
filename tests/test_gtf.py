from pathlib import Path

from ddbj_submission_tools.nuc.ddbj_gff.gtf import parse_gtf

FX = Path(__file__).parent / "fixtures" / "tools"


def test_braker_gtf_builds_hierarchy():
    doc = parse_gtf((FX / "braker3.gtf").read_text())
    g = doc.feature_index["g1"]
    t = doc.feature_index["g1.t1"]
    assert g.type == "gene" and t.type == "transcript" and t.parent_ids == ["g1"]
    assert [c.type for c in t.children].count("CDS") == 1
    cds = next(c for c in t.children if c.type == "CDS")
    assert cds.id == "g1.t1.cds" and all(s.phase is not None for s in cds.spans)
    assert any(c.type == "stop_codon" for c in t.children)
    assert doc.roots == [g]


def test_genemark_gtf_keeps_mrna_type():
    doc = parse_gtf((FX / "genemark_et.gtf").read_text())
    t = doc.feature_index["1_t"]
    assert t.type == "mRNA" and t.parent_ids == ["1_g"]
    assert doc.feature_index["1_g"].type == "gene"
    assert not any(d.code == "gtf-synthesized" for d in doc.diagnostics)


BARE = (
    'chr1\tX\tCDS\t10\t30\t.\t+\t0\tgene_id "ga"; transcript_id "ta";\n'
    'chr1\tX\tCDS\t40\t60\t.\t+\t0\tgene_id "ga"; transcript_id "ta";\n'
)


def test_synthesizes_missing_gene_and_transcript():
    doc = parse_gtf(BARE)
    t = doc.feature_index["ta"]
    g = doc.feature_index["ga"]
    assert t.type == "transcript" and t.parent_ids == ["ga"] and g.type == "gene"
    assert (t.spans[0].start, t.spans[0].end, t.spans[0].strand) == (10, 60, "+")
    assert (g.spans[0].start, g.spans[0].end) == (10, 60)
    assert len(t.children) == 1 and len(t.children[0].spans) == 2
    assert any(d.code == "gtf-synthesized" for d in doc.diagnostics)


def test_stringtie_gtf_synthesizes_gene_and_keeps_attrs():
    doc = parse_gtf((FX / "stringtie.gtf").read_text())
    t = doc.feature_index["STRG.1.1"]
    assert doc.feature_index["STRG.1"].type == "gene"
    assert t.attributes["TPM"] == ["68.092141"]
    assert [c.type for c in t.children] == ["exon"] * len(t.children)
    assert [c.id for c in t.children][:2] == ["STRG.1.1.exon1", "STRG.1.1.exon2"]
