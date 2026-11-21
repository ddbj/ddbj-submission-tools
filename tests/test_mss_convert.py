from Bio.Seq import Seq
from ddbj_submission_tools.nuc.ddbj_gff import parse
from ddbj_submission_tools.nuc.gff2mss.config import MssConfig
from ddbj_submission_tools.nuc.gff2mss.convert import build_entry_features

GFF = """\
##gff-version 3
chr1\tS\tgene\t1\t9\t.\t+\t.\tID=g1;gene=MpX
chr1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g1.1;Parent=g1
chr1\tS\texon\t1\t9\t.\t+\t.\tID=e1;Parent=g1.1
chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g1.1
"""


def cfg():
    return MssConfig(locus_tag_prefix="PFX")


def run(doc, seqs, c):
    diags: list = []
    return build_entry_features(doc, seqs, c, diags), diags


def test_builds_gene_features_per_entry():
    doc = parse(GFF)
    per, diags = run(doc, {"chr1": Seq("ATGAAATAA")}, cfg())
    assert list(per) == ["chr1"]
    keys = [f.key for f in per["chr1"]]
    assert keys == ["mRNA", "CDS"]
    cds = per["chr1"][1]
    assert {q.key: q.value for q in cds.qualifiers}["locus_tag"] == "PFX_000010"


def test_missing_sequence_is_reported_and_skips_entry():
    per, diags = run(parse(GFF), {}, cfg())
    assert any(d.code == "missing-sequence" for d in diags)
    assert per.get("chr1", []) == []

def test_multi_transcript_warns_and_keeps_first():
    gff = GFF + (
        "chr1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g1.2;Parent=g1\n"
        "chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=c2;Parent=g1.2\n"
    )
    doc = parse(gff)
    seqs = {"chr1": Seq("ATGAAATAA")}
    c = cfg(); c.transcript_mode = "minimal"
    per, diags = run(doc, seqs, c)
    assert any(d.code == "multi-transcript" for d in diags)
    # only one mRNA + one CDS kept
    keys = [f.key for f in per["chr1"]]
    assert keys.count("mRNA") == 1 and keys.count("CDS") == 1


def test_default_mode_is_nonredundant_keeps_all_transcripts():
    gff = (
        "##gff-version 3\n"
        "chr1\tS\tgene\t1\t40\t.\t+\t.\tID=g\n"
        "chr1\tS\tmRNA\t1\t30\t.\t+\t.\tID=g.1;Parent=g\n"
        "chr1\tS\texon\t1\t30\t.\t+\t.\tID=e1;Parent=g.1\n"
        "chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g.1\n"
        "chr1\tS\tmRNA\t1\t40\t.\t+\t.\tID=g.2;Parent=g\n"
        "chr1\tS\texon\t1\t40\t.\t+\t.\tID=e2;Parent=g.2\n"
        "chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=c2;Parent=g.2\n"
    )
    doc = parse(gff)
    per, diags = run(doc, {"chr1": Seq("ATGAAATAA" + "C" * 31)}, cfg())
    keys = [f.key for f in per["chr1"]]
    assert keys == ["mRNA", "mRNA", "CDS"]  # nonredundant: 2 mRNA, shared CDS once
    assert not any(d.code == "multi-transcript" for d in diags)  # no warning in nonredundant


def test_locus_tag_unique_across_entries():
    # two seqids, neither gene has a locus_tag attr -> counter must NOT restart per entry
    gff = (
        "##gff-version 3\n"
        "chr1\tS\tgene\t1\t9\t.\t+\t.\tID=a\n"
        "chr1\tS\tmRNA\t1\t9\t.\t+\t.\tID=a.1;Parent=a\n"
        "chr1\tS\texon\t1\t9\t.\t+\t.\tID=ae;Parent=a.1\n"
        "chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=ac;Parent=a.1\n"
        "chr2\tS\tgene\t1\t9\t.\t+\t.\tID=b\n"
        "chr2\tS\tmRNA\t1\t9\t.\t+\t.\tID=b.1;Parent=b\n"
        "chr2\tS\texon\t1\t9\t.\t+\t.\tID=be;Parent=b.1\n"
        "chr2\tS\tCDS\t1\t9\t.\t+\t0\tID=bc;Parent=b.1\n"
    )
    doc = parse(gff)
    seqs = {"chr1": Seq("ATGAAATAA"), "chr2": Seq("ATGAAATAA")}
    per, diags = run(doc, seqs, cfg())
    tags = [q.value for feats in per.values() for f in feats if f.key == "CDS"
            for q in f.qualifiers if q.key == "locus_tag"]
    assert tags == ["PFX_000010", "PFX_000020"]  # monotonic across entries, no duplicate


def test_no_exon_no_cds_mrna_skipped():
    gff = (
        "##gff-version 3\n"
        "chr1\tS\tgene\t1\t9\t.\t+\t.\tID=g\n"
        "chr1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g.1;Parent=g\n"  # no exon/CDS children
    )
    doc = parse(gff)
    per, diags = run(doc, {"chr1": Seq("ATGAAATAA")}, cfg())
    assert any(d.code == "no-exon" for d in diags)
    # no crash; the gene produced no mRNA/CDS features
    assert not any(f.key in ("mRNA", "CDS") for feats in per.values() for f in feats)
