from ddbj_submission_tools.nuc.ddbj_gff import parse
from ddbj_submission_tools.nuc.ddbj_gff.repair import REGISTRY, RepairContext, run_detect, run_apply, apply_all_names
from ddbj_submission_tools.nuc.ddbj_gff.repair.operations import _duplicate_mrna_losers, _is_five_prime_truncation

HDR = "##gff-version 3\n##sequence-region c 1 100000\n"


def _doc(body):
    return parse(HDR + body)


# --- plain 5' truncation, plus strand: t2 starts later, same 3' end -------
DUP_PLUS = (
    "c\tx\tgene\t100\t600\t.\t+\t.\tID=g\n"
    "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.t1;Parent=g\n"
    "c\tx\texon\t100\t200\t.\t+\t.\tID=g.t1.e1;Parent=g.t1\n"
    "c\tx\texon\t500\t600\t.\t+\t.\tID=g.t1.e2;Parent=g.t1\n"
    "c\tx\tCDS\t120\t200\t.\t+\t0\tID=g.t1.c1;Parent=g.t1\n"
    "c\tx\tCDS\t500\t580\t.\t+\t0\tID=g.t1.c2;Parent=g.t1\n"
    "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.t2;Parent=g\n"
    "c\tx\texon\t100\t200\t.\t+\t.\tID=g.t2.e1;Parent=g.t2\n"
    "c\tx\texon\t500\t600\t.\t+\t.\tID=g.t2.e2;Parent=g.t2\n"
    "c\tx\tCDS\t150\t200\t.\t+\t0\tID=g.t2.c1;Parent=g.t2\n"
    "c\tx\tCDS\t500\t580\t.\t+\t0\tID=g.t2.c2;Parent=g.t2\n"
)


def test_truncation_helper_accepts_identical_cds():
    segs = [(120, 200), (500, 580)]
    assert _is_five_prime_truncation(segs, segs, "+") is True


def test_truncation_helper_rejects_new_splice_junction():
    # S is contained in L but introduces an intron L does not have
    assert _is_five_prime_truncation([(120, 200), (300, 500)], [(120, 500)], "+") is False


def test_truncation_helper_rejects_three_prime_difference():
    # same start, shorter 3' end on the plus strand -> not a 5' truncation
    assert _is_five_prime_truncation([(120, 400)], [(120, 500)], "+") is False


def test_losers_plain_five_prime_truncation():
    losers = _duplicate_mrna_losers(_doc(DUP_PLUS))
    assert losers == {"g.t2": ("g.t1", 132, 162)}


def test_losers_g200_shape_dropped_segment():
    # t2's start codon moves into a later exon, so its CDS loses a whole segment
    doc = _doc(
        "c\tx\tgene\t100\t900\t.\t+\t.\tID=h\n"
        "c\tx\tmRNA\t100\t900\t.\t+\t.\tID=h.t1;Parent=h\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=h.t1.e1;Parent=h.t1\n"
        "c\tx\texon\t300\t400\t.\t+\t.\tID=h.t1.e2;Parent=h.t1\n"
        "c\tx\texon\t700\t900\t.\t+\t.\tID=h.t1.e3;Parent=h.t1\n"
        "c\tx\tCDS\t120\t200\t.\t+\t0\tID=h.t1.c1;Parent=h.t1\n"
        "c\tx\tCDS\t300\t400\t.\t+\t0\tID=h.t1.c2;Parent=h.t1\n"
        "c\tx\tCDS\t700\t800\t.\t+\t0\tID=h.t1.c3;Parent=h.t1\n"
        "c\tx\tmRNA\t100\t900\t.\t+\t.\tID=h.t2;Parent=h\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=h.t2.e1;Parent=h.t2\n"
        "c\tx\texon\t300\t400\t.\t+\t.\tID=h.t2.e2;Parent=h.t2\n"
        "c\tx\texon\t700\t900\t.\t+\t.\tID=h.t2.e3;Parent=h.t2\n"
        "c\tx\tCDS\t350\t400\t.\t+\t0\tID=h.t2.c1;Parent=h.t2\n"
        "c\tx\tCDS\t700\t800\t.\t+\t0\tID=h.t2.c2;Parent=h.t2\n"
    )
    assert _duplicate_mrna_losers(doc) == {"h.t2": ("h.t1", 152, 283)}


def test_losers_minus_strand_truncates_at_genomic_high_end():
    # on '-' the 5' end is the genomic high end; the 3' terminus is the genomic start
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t-\t.\tID=k\n"
        "c\tx\tmRNA\t100\t600\t.\t-\t.\tID=k.t1;Parent=k\n"
        "c\tx\texon\t100\t200\t.\t-\t.\tID=k.t1.e1;Parent=k.t1\n"
        "c\tx\texon\t500\t600\t.\t-\t.\tID=k.t1.e2;Parent=k.t1\n"
        "c\tx\tCDS\t120\t200\t.\t-\t0\tID=k.t1.c1;Parent=k.t1\n"
        "c\tx\tCDS\t500\t580\t.\t-\t0\tID=k.t1.c2;Parent=k.t1\n"
        "c\tx\tmRNA\t100\t600\t.\t-\t.\tID=k.t2;Parent=k\n"
        "c\tx\texon\t100\t200\t.\t-\t.\tID=k.t2.e1;Parent=k.t2\n"
        "c\tx\texon\t500\t600\t.\t-\t.\tID=k.t2.e2;Parent=k.t2\n"
        "c\tx\tCDS\t120\t200\t.\t-\t0\tID=k.t2.c1;Parent=k.t2\n"
        "c\tx\tCDS\t500\t550\t.\t-\t0\tID=k.t2.c2;Parent=k.t2\n"
    )
    assert _duplicate_mrna_losers(doc) == {"k.t2": ("k.t1", 132, 162)}


def test_losers_identical_cds_tie_breaks_on_id():
    doc = _doc(DUP_PLUS.replace(
        "c\tx\tCDS\t150\t200\t.\t+\t0\tID=g.t2.c1;Parent=g.t2\n",
        "c\tx\tCDS\t120\t200\t.\t+\t0\tID=g.t2.c1;Parent=g.t2\n"))
    assert _duplicate_mrna_losers(doc) == {"g.t2": ("g.t1", 162, 162)}


def test_losers_group_of_three():
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=m\n"
        + "".join(
            f"c\tx\tmRNA\t100\t600\t.\t+\t.\tID=m.t{i};Parent=m\n"
            f"c\tx\texon\t100\t600\t.\t+\t.\tID=m.t{i}.e1;Parent=m.t{i}\n"
            f"c\tx\tCDS\t{start}\t580\t.\t+\t0\tID=m.t{i}.c1;Parent=m.t{i}\n"
            for i, start in ((1, 120), (2, 150), (3, 180))))
    assert _duplicate_mrna_losers(doc) == {
        "m.t2": ("m.t1", 431, 461),
        "m.t3": ("m.t1", 401, 461),
    }


def test_losers_divergent_group_untouched():
    # t2 adds a splice junction t1 does not have -> divergent
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=d\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=d.t1;Parent=d\n"
        "c\tx\texon\t100\t600\t.\t+\t.\tID=d.t1.e1;Parent=d.t1\n"
        "c\tx\tCDS\t120\t500\t.\t+\t0\tID=d.t1.c1;Parent=d.t1\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=d.t2;Parent=d\n"
        "c\tx\texon\t100\t600\t.\t+\t.\tID=d.t2.e1;Parent=d.t2\n"
        "c\tx\tCDS\t120\t200\t.\t+\t0\tID=d.t2.c1;Parent=d.t2\n"
        "c\tx\tCDS\t300\t500\t.\t+\t0\tID=d.t2.c2;Parent=d.t2\n"
    )
    assert _duplicate_mrna_losers(doc) == {}


def test_losers_skips_group_with_a_cds_less_member():
    doc = _doc(DUP_PLUS.replace(
        "c\tx\tCDS\t150\t200\t.\t+\t0\tID=g.t2.c1;Parent=g.t2\n"
        "c\tx\tCDS\t500\t580\t.\t+\t0\tID=g.t2.c2;Parent=g.t2\n", ""))
    assert _duplicate_mrna_losers(doc) == {}


def test_losers_skips_trans_spliced():
    doc = _doc(DUP_PLUS.replace(
        "ID=g.t2;Parent=g\n", "ID=g.t2;Parent=g;exception=trans-splicing\n"))
    assert _duplicate_mrna_losers(doc) == {}


def test_losers_does_not_group_across_genes():
    # same exon structure, but two different genes -> not a duplicate group
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=p\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=p.t1;Parent=p\n"
        "c\tx\texon\t100\t600\t.\t+\t.\tID=p.t1.e1;Parent=p.t1\n"
        "c\tx\tCDS\t120\t580\t.\t+\t0\tID=p.t1.c1;Parent=p.t1\n"
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=q\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=q.t1;Parent=q\n"
        "c\tx\texon\t100\t600\t.\t+\t.\tID=q.t1.e1;Parent=q.t1\n"
        "c\tx\tCDS\t150\t580\t.\t+\t0\tID=q.t1.c1;Parent=q.t1\n"
    )
    assert _duplicate_mrna_losers(doc) == {}


DUP = "duplicate-mrna-keep-longest-cds"


def test_detect_reports_the_shorter_mrna():
    cands = run_detect(_doc(DUP_PLUS), RepairContext(), [DUP])
    assert len(cands) == 1
    c = cands[0]
    assert c.operation == DUP and c.feature_id == "g.t2"
    assert c.payload == {"keep": "g.t1", "cds_len": 132, "keep_cds_len": 162}
    assert "132 bp vs 162 bp" in c.detail


def test_apply_removes_the_shorter_mrna_and_its_children():
    doc = _doc(DUP_PLUS)
    changes = run_apply(doc, RepairContext(), [DUP])
    assert doc.feature_index.get("g.t2") is None
    for child in ("g.t2.e1", "g.t2.e2", "g.t2.c1", "g.t2.c2"):
        assert doc.feature_index.get(child) is None
    keep = doc.feature_index.get("g.t1")
    assert keep is not None                                  # winner untouched
    assert {c.id for c in keep.children} == {"g.t1.e1", "g.t1.e2", "g.t1.c1", "g.t1.c2"}
    gene = doc.feature_index.get("g")
    assert [c.id for c in gene.children] == ["g.t1"]         # gene keeps the winner
    assert [c.action for c in changes] == ["drop-duplicate-mrna"]
    assert "kept 'g.t1'" in changes[0].message


def test_apply_is_idempotent():
    doc = _doc(DUP_PLUS)
    run_apply(doc, RepairContext(), [DUP])
    assert run_detect(doc, RepairContext(), [DUP]) == []


def test_apply_leaves_divergent_group_alone():
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=d\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=d.t1;Parent=d\n"
        "c\tx\texon\t100\t600\t.\t+\t.\tID=d.t1.e1;Parent=d.t1\n"
        "c\tx\tCDS\t120\t500\t.\t+\t0\tID=d.t1.c1;Parent=d.t1\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=d.t2;Parent=d\n"
        "c\tx\texon\t100\t600\t.\t+\t.\tID=d.t2.e1;Parent=d.t2\n"
        "c\tx\tCDS\t120\t200\t.\t+\t0\tID=d.t2.c1;Parent=d.t2\n"
        "c\tx\tCDS\t300\t500\t.\t+\t0\tID=d.t2.c2;Parent=d.t2\n"
    )
    assert run_apply(doc, RepairContext(), [DUP]) == []
    assert doc.feature_index.get("d.t2") is not None


def test_apply_keeps_child_shared_with_the_surviving_mrna():
    # the shared exon belongs to both transcripts; it must survive with g.t2 detached
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=g\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.t1;Parent=g\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.t2;Parent=g\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=g.shared;Parent=g.t1,g.t2\n"
        "c\tx\texon\t500\t600\t.\t+\t.\tID=g.t1.e2;Parent=g.t1\n"
        "c\tx\texon\t500\t600\t.\t+\t.\tID=g.t2.e2;Parent=g.t2\n"
        "c\tx\tCDS\t120\t200\t.\t+\t0\tID=g.t1.c1;Parent=g.t1\n"
        "c\tx\tCDS\t500\t580\t.\t+\t0\tID=g.t1.c2;Parent=g.t1\n"
        "c\tx\tCDS\t150\t200\t.\t+\t0\tID=g.t2.c1;Parent=g.t2\n"
        "c\tx\tCDS\t500\t580\t.\t+\t0\tID=g.t2.c2;Parent=g.t2\n"
    )
    changes = run_apply(doc, RepairContext(), [DUP])
    assert doc.feature_index.get("g.t2") is None
    shared = doc.feature_index.get("g.shared")
    assert shared is not None and shared.parent_ids == ["g.t1"]
    assert shared.attributes["Parent"] == ["g.t1"]
    assert any("shared with another parent" in c.message for c in changes)


def test_operation_is_destructive_and_excluded_from_apply_all():
    assert REGISTRY[DUP].destructive is True
    assert REGISTRY[DUP].requires_sequence is False
    assert DUP not in apply_all_names()


def test_apply_group_of_three_keeps_only_the_longest():
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=m\n"
        + "".join(
            f"c\tx\tmRNA\t100\t600\t.\t+\t.\tID=m.t{i};Parent=m\n"
            f"c\tx\texon\t100\t600\t.\t+\t.\tID=m.t{i}.e1;Parent=m.t{i}\n"
            f"c\tx\tCDS\t{start}\t580\t.\t+\t0\tID=m.t{i}.c1;Parent=m.t{i}\n"
            for i, start in ((1, 120), (2, 150), (3, 180))))
    changes = run_apply(doc, RepairContext(), [DUP])
    assert doc.feature_index.get("m.t1") is not None
    assert doc.feature_index.get("m.t2") is None
    assert doc.feature_index.get("m.t3") is None
    assert [c.id for c in doc.feature_index["m"].children] == ["m.t1"]
    assert len(changes) == 2


def test_apply_minus_strand_keeps_the_longer_cds():
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t-\t.\tID=k\n"
        "c\tx\tmRNA\t100\t600\t.\t-\t.\tID=k.t1;Parent=k\n"
        "c\tx\texon\t100\t200\t.\t-\t.\tID=k.t1.e1;Parent=k.t1\n"
        "c\tx\texon\t500\t600\t.\t-\t.\tID=k.t1.e2;Parent=k.t1\n"
        "c\tx\tCDS\t120\t200\t.\t-\t0\tID=k.t1.c1;Parent=k.t1\n"
        "c\tx\tCDS\t500\t580\t.\t-\t0\tID=k.t1.c2;Parent=k.t1\n"
        "c\tx\tmRNA\t100\t600\t.\t-\t.\tID=k.t2;Parent=k\n"
        "c\tx\texon\t100\t200\t.\t-\t.\tID=k.t2.e1;Parent=k.t2\n"
        "c\tx\texon\t500\t600\t.\t-\t.\tID=k.t2.e2;Parent=k.t2\n"
        "c\tx\tCDS\t120\t200\t.\t-\t0\tID=k.t2.c1;Parent=k.t2\n"
        "c\tx\tCDS\t500\t550\t.\t-\t0\tID=k.t2.c2;Parent=k.t2\n"
    )
    run_apply(doc, RepairContext(), [DUP])
    assert doc.feature_index.get("k.t1") is not None      # 162 bp CDS kept
    assert doc.feature_index.get("k.t2") is None          # 132 bp CDS dropped


def test_apply_twice_is_a_no_op():
    doc = _doc(DUP_PLUS)
    assert len(run_apply(doc, RepairContext(), [DUP])) == 1
    assert run_apply(doc, RepairContext(), [DUP]) == []


def test_strandless_mrna_uses_the_forward_rule():
    # strand '.' must follow the forward convention: a shorter genomic END is a 3'
    # truncation (a different product), not a 5' one, so nothing may be collapsed.
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t.\t.\tID=s\n"
        "c\tx\tmRNA\t100\t600\t.\t.\t.\tID=s.t1;Parent=s\n"
        "c\tx\texon\t100\t600\t.\t.\t.\tID=s.t1.e1;Parent=s.t1\n"
        "c\tx\tCDS\t120\t580\t.\t.\t0\tID=s.t1.c1;Parent=s.t1\n"
        "c\tx\tmRNA\t100\t600\t.\t.\t.\tID=s.t2;Parent=s\n"
        "c\tx\texon\t100\t600\t.\t.\t.\tID=s.t2.e1;Parent=s.t2\n"
        "c\tx\tCDS\t120\t400\t.\t.\t0\tID=s.t2.c1;Parent=s.t2\n"
    )
    assert run_detect(doc, RepairContext(), [DUP]) == []
    assert run_apply(doc, RepairContext(), [DUP]) == []
    assert doc.feature_index.get("s.t2") is not None
