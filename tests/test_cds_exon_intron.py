from ddbj_submission_tools.nuc.ddbj_gff import parse
from ddbj_submission_tools.nuc.ddbj_gff.model import cds_exon_intron_conflict

HDR = "##gff-version 3\n##sequence-region c 1 100000\n"


def _mrna(doc, mid):
    return doc.feature_index[mid]


def _g(rows):
    return HDR + "".join(rows)


def test_conflict_interior_utr_exon():
    # exons 100-200,300-400,500-600 ; CDS 100-200,500-600 (skips the middle exon)
    gff = _g([
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=g\n",
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.m;Parent=g\n",
        "c\tx\texon\t100\t200\t.\t+\t.\tID=g.e1;Parent=g.m\n",
        "c\tx\texon\t300\t400\t.\t+\t.\tID=g.e2;Parent=g.m\n",
        "c\tx\texon\t500\t600\t.\t+\t.\tID=g.e3;Parent=g.m\n",
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g.c1;Parent=g.m\n",
        "c\tx\tCDS\t500\t600\t.\t+\t0\tID=g.c2;Parent=g.m\n",
    ])
    conflicts = cds_exon_intron_conflict(_mrna(parse(gff), "g.m"))
    assert conflicts == [(201, 499)]


def test_clean_consistent_gene():
    # CDS within exons, same junctions; terminal 5'/3' UTR
    gff = _g([
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=g;\n",
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.m;Parent=g\n",
        "c\tx\texon\t100\t200\t.\t+\t.\tID=g.e1;Parent=g.m\n",
        "c\tx\texon\t300\t400\t.\t+\t.\tID=g.e2;Parent=g.m\n",
        "c\tx\texon\t500\t600\t.\t+\t.\tID=g.e3;Parent=g.m\n",
        "c\tx\tCDS\t150\t200\t.\t+\t0\tID=g.c1;Parent=g.m\n",
        "c\tx\tCDS\t300\t400\t.\t+\t0\tID=g.c2;Parent=g.m\n",
        "c\tx\tCDS\t500\t550\t.\t+\t0\tID=g.c3;Parent=g.m\n",
    ])
    assert cds_exon_intron_conflict(_mrna(parse(gff), "g.m")) == []


def test_terminal_multiexon_utr_not_flagged():
    # 5'UTR-only exon 50-150, then coding exons; 3'UTR-only exon 900-1000
    gff = _g([
        "c\tx\tgene\t50\t1000\t.\t+\t.\tID=g\n",
        "c\tx\tmRNA\t50\t1000\t.\t+\t.\tID=g.m;Parent=g\n",
        "c\tx\texon\t50\t150\t.\t+\t.\tID=g.e0;Parent=g.m\n",
        "c\tx\texon\t300\t400\t.\t+\t.\tID=g.e1;Parent=g.m\n",
        "c\tx\texon\t600\t700\t.\t+\t.\tID=g.e2;Parent=g.m\n",
        "c\tx\texon\t900\t1000\t.\t+\t.\tID=g.e3;Parent=g.m\n",
        "c\tx\tCDS\t300\t400\t.\t+\t0\tID=g.c1;Parent=g.m\n",
        "c\tx\tCDS\t600\t700\t.\t+\t0\tID=g.c2;Parent=g.m\n",
    ])
    assert cds_exon_intron_conflict(_mrna(parse(gff), "g.m")) == []


def test_single_cds_no_conflict():
    gff = _g([
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=g\n",
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.m;Parent=g\n",
        "c\tx\texon\t100\t600\t.\t+\t.\tID=g.e1;Parent=g.m\n",
        "c\tx\tCDS\t100\t600\t.\t+\t0\tID=g.c1;Parent=g.m\n",
    ])
    assert cds_exon_intron_conflict(_mrna(parse(gff), "g.m")) == []


from ddbj_submission_tools.nuc.ddbj_gff.validate import validate


def test_rule_flags_internal_utr():
    gff = _g([
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=g\n",
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.m;Parent=g\n",
        "c\tx\texon\t100\t200\t.\t+\t.\tID=g.e1;Parent=g.m\n",
        "c\tx\texon\t300\t400\t.\t+\t.\tID=g.e2;Parent=g.m\n",
        "c\tx\texon\t500\t600\t.\t+\t.\tID=g.e3;Parent=g.m\n",
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g.c1;Parent=g.m\n",
        "c\tx\tCDS\t500\t600\t.\t+\t0\tID=g.c2;Parent=g.m\n",
    ])
    diags = [d for d in validate(parse(gff)) if d.code == "cds-exon-intron-mismatch"]
    assert len(diags) == 1
    assert diags[0].severity.value == "WARNING"


def test_rule_silent_on_clean_gene():
    gff = _g([
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=g\n",
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.m;Parent=g\n",
        "c\tx\texon\t100\t200\t.\t+\t.\tID=g.e1;Parent=g.m\n",
        "c\tx\texon\t500\t600\t.\t+\t.\tID=g.e2;Parent=g.m\n",
        "c\tx\tCDS\t150\t200\t.\t+\t0\tID=g.c1;Parent=g.m\n",
        "c\tx\tCDS\t500\t550\t.\t+\t0\tID=g.c2;Parent=g.m\n",
    ])
    assert not [d for d in validate(parse(gff)) if d.code == "cds-exon-intron-mismatch"]
