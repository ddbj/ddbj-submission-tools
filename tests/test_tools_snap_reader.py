from pathlib import Path

from ddbj_submission_tools.nuc.ddbj_gff.tools.snap_reader import parse_snap

FX = Path(__file__).parent / "fixtures" / "tools"


def test_groups_rows_by_name_into_gene_mrna_cds():
    text = (FX / "snap.gff3").read_text()
    rows = [l for l in text.splitlines() if l and not l.startswith("#")]
    names = []
    for r in rows:
        n = r.split("\t")[7].split("=", 1)[1]
        if n not in names:
            names.append(n)
    doc = parse_snap(text)
    genes = [f for f in doc.roots]
    assert [g.id for g in genes] == names and all(g.type == "gene" for g in genes)
    g = doc.feature_index[names[0]]
    m = doc.feature_index[names[0] + ".t1"]
    assert m.type == "mRNA" and m.parent_ids == [names[0]] and g.children == [m]
    cds = doc.feature_index[names[0] + ".t1.cds"]
    assert cds.type == "CDS" and cds.parent_ids == [names[0] + ".t1"]
    assert len(cds.spans) == sum(1 for r in rows if r.endswith("Name=" + names[0]))
    assert all(s.phase is None for s in cds.spans)
    assert cds.spans[0].strand == rows[0].split("\t")[6]
    assert (g.spans[0].start, g.spans[0].end) == (min(s.start for s in cds.spans), max(s.end for s in cds.spans))


def test_nine_column_rows_are_tolerated():
    text = "chr1\tsnap\tCDS\t10\t20\t.\t+\t.\tName=x.1\nchr1\tsnap\tCDS\t30\t40\t.\t+\t.\tName=x.1\n"
    doc = parse_snap(text)
    assert len(doc.feature_index["x.1.t1.cds"].spans) == 2
