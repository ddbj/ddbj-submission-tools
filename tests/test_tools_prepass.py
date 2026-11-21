from ddbj_submission_tools.nuc.ddbj_gff.tools.prepass import dedupe_ids_per_seqid

GFF = """##gff-version 3
chr1\tA\tgene\t1\t9\t.\t+\t.\tID=g1
chr1\tA\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1
chr1\tA\tCDS\t1\t3\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1
chr1\tA\tCDS\t7\t9\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1
chr2\tA\tgene\t1\t9\t.\t+\t.\tID=g1
chr2\tA\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1
##FASTA
>chr1
ID=g1
"""


def test_renames_cross_seqid_ids_and_parents():
    out, n = dedupe_ids_per_seqid(GFF)
    assert n == 2                     # g1 and g1.t1; g1.t1.cds only on chr1
    assert "ID=chr1.g1\n" in out and "ID=chr2.g1\n" in out
    assert "ID=chr2.g1.t1;Parent=chr2.g1" in out
    assert out.count("ID=g1.t1.cds") == 2          # same-seqid repeat untouched
    assert "Parent=chr1.g1.t1" in out
    assert out.endswith("##FASTA\n>chr1\nID=g1\n")  # FASTA section untouched


def test_no_collision_is_identity():
    text = GFF.split("chr2")[0]
    assert dedupe_ids_per_seqid(text) == (text, 0)


def test_gtf_ids():
    gtf = ('chr1\tA\tgene\t1\t9\t.\t+\t.\tg1\n'
           'chr1\tA\tCDS\t1\t9\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
           'chr2\tA\tgene\t1\t9\t.\t+\t.\tg1\n'
           'chr2\tA\tCDS\t1\t9\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n')
    out, n = dedupe_ids_per_seqid(gtf)
    assert n == 2 and 'gene_id "chr2.g1"' in out and out.splitlines()[2].endswith("\tchr2.g1")
    assert 'transcript_id "chr1.g1.t1"' in out


def test_gtf_detection_ignores_comment_lines():
    gff = ('# made from gene_id "x" style input\n'
           'chr1\tA\tgene\t1\t9\t.\t+\t.\tID=g1\nchr2\tA\tgene\t1\t9\t.\t+\t.\tID=g1\n')
    out, n = dedupe_ids_per_seqid(gff)
    assert n == 1 and "ID=chr1.g1" in out and "chr1.ID=" not in out
