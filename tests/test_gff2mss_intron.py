"""gff2mss writes intron features only on request ([transcript] emit_intron), with the gene's locus_tag."""
from Bio.Seq import Seq

from ddbj_gff import parse
from gff2mss.config import MssConfig, load_config
from tests.mss_helpers import ann_text

GFF = """##gff-version 3
chr1\tA\tgene\t1\t30\t.\t+\t.\tID=g1
chr1\tA\tmRNA\t1\t30\t.\t+\t.\tID=g1.t1;Parent=g1
chr1\tA\texon\t1\t9\t.\t+\t.\tParent=g1.t1
chr1\tA\tintron\t10\t19\t.\t+\t.\tParent=g1.t1
chr1\tA\texon\t20\t30\t.\t+\t.\tParent=g1.t1
chr1\tA\tCDS\t1\t9\t.\t+\t0\tParent=g1.t1
chr1\tA\tCDS\t20\t22\t.\t+\t0\tParent=g1.t1
"""
SEQS = {"chr1": Seq("ATGAAACCC" + "GTAAGTTTAG" + "TAA" + "A" * 8)}


def _ann(**kw):
    return ann_text(parse(GFF), SEQS, MssConfig(locus_tag_prefix="T", **kw))


def _rows(ann, key):
    return [r.split("\t") for r in ann.splitlines() if r.split("\t")[1:2] == [key]]


def test_introns_are_not_written_by_default():
    ann = _ann()
    assert not _rows(ann, "intron") and _rows(ann, "CDS")


def test_emit_intron_writes_them_with_the_gene_locus_tag():
    ann = _ann(emit_intron=True)
    (intron,) = _rows(ann, "intron")
    (cds,) = _rows(ann, "CDS")
    assert intron[2] == "10..19" and intron[3:5] == ["locus_tag", cds[4]]


def test_config_emit_intron(tmp_path):
    p = tmp_path / "m.toml"
    p.write_text('[locus_tag]\nprefix="T"\n')
    assert load_config(str(p))[0].emit_intron is False
    p.write_text('[locus_tag]\nprefix="T"\n[transcript]\nemit_intron = true\n')
    assert load_config(str(p))[0].emit_intron is True
