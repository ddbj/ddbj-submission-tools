from ddbj_gff import parse
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.normalize.normalize import normalize

GFF = """##gff-version 3
chr1\tGlimmerHMM\tmRNA\t2929\t3249\t.\t-\t.\tID=chr1.path1.gene1;Name=chr1.path1.gene1
chr1\tGlimmerHMM\tCDS\t2929\t3000\t.\t-\t0\tID=chr1.cds1.1;Parent=chr1.path1.gene1
chr1\tGlimmerHMM\tCDS\t3100\t3249\t.\t-\t2\tID=chr1.cds1.2;Parent=chr1.path1.gene1
"""


def test_wraps_parentless_mrna_in_gene():
    norm, rep = normalize(parse(GFF), config=NormalizeConfig(taxid=1, wrap_mrna_in_gene=True))
    gene = norm.feature_index["chr1.path1.gene1.gene"]
    mrna = norm.feature_index["chr1.path1.gene1"]
    assert gene.type == "gene"
    assert (gene.spans[0].start, gene.spans[0].end, gene.spans[0].strand) == (2929, 3249, "-")
    assert norm.roots == [gene]
    assert mrna.parent_ids == ["chr1.path1.gene1.gene"] and mrna.parents == [gene]
    assert mrna.attributes["Parent"] == ["chr1.path1.gene1.gene"]
    assert gene.children == [mrna]
    assert [(c.action, c.target) for c in rep.applied if c.action == "add-parent-feature"] == \
        [("add-parent-feature", "chr1.path1.gene1.gene")]


def test_default_off_keeps_mrna_as_root():
    norm, _ = normalize(parse(GFF), config=NormalizeConfig(taxid=1))
    assert [r.type for r in norm.roots] == ["mRNA"]


def test_wrap_is_linear_in_root_count():
    import time
    rows = "".join(f"chr1\tG\tmRNA\t{i*10}\t{i*10+5}\t.\t+\t.\tID=m{i}\n"
                   f"chr1\tG\tCDS\t{i*10}\t{i*10+5}\t.\t+\t0\tID=m{i}.c;Parent=m{i}\n" for i in range(20000))
    doc = parse("##gff-version 3\n" + rows)
    t0 = time.time()
    norm, _ = normalize(doc, config=NormalizeConfig(taxid=1, wrap_mrna_in_gene=True))
    assert len(norm.roots) == 20000 and time.time() - t0 < 5
