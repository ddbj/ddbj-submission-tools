from ddbj_gff import parse
from ddbj_gff.repair import REGISTRY, RepairContext, run_detect, run_apply, apply_all_names

HDR = "##gff-version 3\n##sequence-region c 1 100000\n"
OP = "short-intron-artificial-location"
VALUE = "low-quality sequence region"


def _doc(body):
    return parse(HDR + body)


def _al(doc, fid):
    f = doc.feature_index.get(fid)
    return None if f is None else f.attributes.get("artificial_location")


# exon-side 8 bp intron (the g6682.t1 shape): exons 100-200 and 209-600, CDS contiguous
EXON_SHORT = (
    "c\tx\tgene\t100\t600\t.\t+\t.\tID=g\n"
    "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.m;Parent=g\n"
    "c\tx\texon\t100\t200\t.\t+\t.\tID=g.e1;Parent=g.m\n"
    "c\tx\texon\t209\t600\t.\t+\t.\tID=g.e2;Parent=g.m\n"
    "c\tx\tCDS\t120\t200\t.\t+\t0\tID=g.c1;Parent=g.m\n"
    "c\tx\tCDS\t209\t580\t.\t+\t0\tID=g.c2;Parent=g.m\n"
)


def test_detect_marks_the_mrna_for_an_exon_side_short_intron():
    cands = run_detect(_doc(EXON_SHORT), RepairContext(), [OP])
    ids = sorted(c.feature_id for c in cands)
    assert ids == ["g.c1", "g.c2", "g.m"]      # CDS gap is also 8 bp here
    assert all(c.payload == {"shortest_intron": 8, "value": VALUE} for c in cands)
    assert "8 bp" in cands[0].detail


def test_detect_mrna_only_when_cds_is_unaffected():
    # exon gap 201..208 = 8 bp; the CDS stops before it, so CDS has no intron at all
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=h\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=h.m;Parent=h\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=h.e1;Parent=h.m\n"
        "c\tx\texon\t209\t600\t.\t+\t.\tID=h.e2;Parent=h.m\n"
        "c\tx\tCDS\t120\t180\t.\t+\t0\tID=h.c1;Parent=h.m\n"
    )
    assert [c.feature_id for c in run_detect(doc, RepairContext(), [OP])] == ["h.m"]


def test_detect_marks_every_cds_child_for_a_cds_side_short_intron():
    # exons contiguous (no exon intron); the CDS has a 5 bp gap 201..205
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=k\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=k.m;Parent=k\n"
        "c\tx\texon\t100\t600\t.\t+\t.\tID=k.e1;Parent=k.m\n"
        "c\tx\tCDS\t120\t200\t.\t+\t0\tID=k.c1;Parent=k.m\n"
        "c\tx\tCDS\t206\t580\t.\t+\t0\tID=k.c2;Parent=k.m\n"
    )
    assert sorted(c.feature_id for c in run_detect(doc, RepairContext(), [OP])) == ["k.c1", "k.c2"]


def test_detect_uses_an_explicit_short_intron_child():
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=n\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=n.m;Parent=n\n"
        "c\tx\texon\t100\t600\t.\t+\t.\tID=n.e1;Parent=n.m\n"
        "c\tx\tintron\t300\t302\t.\t+\t.\tID=n.i1;Parent=n.m\n"
        "c\tx\tCDS\t120\t580\t.\t+\t0\tID=n.c1;Parent=n.m\n"
    )
    cands = run_detect(doc, RepairContext(), [OP])
    assert [c.feature_id for c in cands] == ["n.m"]
    assert cands[0].payload["shortest_intron"] == 3


def test_detect_silent_when_all_introns_are_long_enough():
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=p\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=p.m;Parent=p\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=p.e1;Parent=p.m\n"
        "c\tx\texon\t211\t600\t.\t+\t.\tID=p.e2;Parent=p.m\n"   # gap 201..210 = 10 bp
        "c\tx\tCDS\t120\t200\t.\t+\t0\tID=p.c1;Parent=p.m\n"
        "c\tx\tCDS\t211\t580\t.\t+\t0\tID=p.c2;Parent=p.m\n"
    )
    assert run_detect(doc, RepairContext(), [OP]) == []


def test_detect_skips_features_that_already_declare_an_exception():
    doc = _doc(EXON_SHORT
               .replace("ID=g.m;Parent=g\n", "ID=g.m;Parent=g;artificial_location=low-quality sequence region\n")
               .replace("ID=g.c1;Parent=g.m\n", "ID=g.c1;Parent=g.m;ribosomal_slippage=true\n")
               .replace("ID=g.c2;Parent=g.m\n", "ID=g.c2;Parent=g.m;artificial_location=heterogeneous population sequenced\n"))
    assert run_detect(doc, RepairContext(), [OP]) == []


def test_detect_skips_trans_spliced_mrna():
    doc = _doc(EXON_SHORT.replace("ID=g.m;Parent=g\n", "ID=g.m;Parent=g;exception=trans-splicing\n"))
    assert run_detect(doc, RepairContext(), [OP]) == []


def test_apply_writes_the_fixed_value():
    doc = _doc(EXON_SHORT)
    changes = run_apply(doc, RepairContext(), [OP])
    assert _al(doc, "g.m") == [VALUE]
    assert _al(doc, "g.c1") == [VALUE] and _al(doc, "g.c2") == [VALUE]
    assert {c.action for c in changes} == {"add-artificial-location"}


def test_apply_does_not_change_the_feature_graph():
    doc = _doc(EXON_SHORT)
    before = ([f.id for f in doc.features], sorted(doc.feature_index),
              [r.id for r in doc.roots],
              {f.id: sorted(c.id for c in f.children) for f in doc.features})
    run_apply(doc, RepairContext(), [OP])
    after = ([f.id for f in doc.features], sorted(doc.feature_index),
             [r.id for r in doc.roots],
             {f.id: sorted(c.id for c in f.children) for f in doc.features})
    assert before == after


def test_apply_is_idempotent():
    doc = _doc(EXON_SHORT)
    run_apply(doc, RepairContext(), [OP])
    assert run_detect(doc, RepairContext(), [OP]) == []
    assert run_apply(doc, RepairContext(), [OP]) == []


def test_operation_is_destructive_and_excluded_from_apply_all():
    assert REGISTRY[OP].destructive is True
    assert REGISTRY[OP].requires_sequence is False
    assert OP not in apply_all_names()


def test_exonless_mrna_is_marked_from_its_cds_spans():
    # pass_wrap_cds_in_mrna produces mRNAs with no exon children; gff2mss then builds the
    # mRNA location from the CDS spans, so a short CDS gap shows up on the .ann mRNA too.
    doc = _doc(
        "c\tx\tgene\t2000\t2600\t.\t+\t.\tID=w\n"
        "c\tx\tmRNA\t2000\t2600\t.\t+\t.\tID=w.m;Parent=w\n"
        "c\tx\tCDS\t2000\t2101\t.\t+\t0\tID=w.c1;Parent=w.m\n"
        "c\tx\tCDS\t2109\t2600\t.\t+\t0\tID=w.c2;Parent=w.m\n"
    )
    assert sorted(c.feature_id for c in run_detect(doc, RepairContext(), [OP])) == [
        "w.c1", "w.c2", "w.m"]
    run_apply(doc, RepairContext(), [OP])
    assert _al(doc, "w.m") == [VALUE]
    assert _al(doc, "w.c1") == [VALUE] and _al(doc, "w.c2") == [VALUE]


def test_skips_when_only_a_cds_child_is_trans_spliced():
    # normalize marks the CDS, not the mRNA, so the guard must look at the children too
    doc = _doc(
        "c\tx\tgene\t2000\t2600\t.\t+\t.\tID=t\n"
        "c\tx\tmRNA\t2000\t2600\t.\t+\t.\tID=t.m;Parent=t\n"
        "c\tx\texon\t2000\t2101\t.\t+\t.\tID=t.e1;Parent=t.m\n"
        "c\tx\texon\t2109\t2600\t.\t+\t.\tID=t.e2;Parent=t.m\n"
        "c\tx\tCDS\t2000\t2101\t.\t+\t0\tID=t.c1;Parent=t.m;exception=trans-splicing\n"
        "c\tx\tCDS\t2109\t2600\t.\t+\t0\tID=t.c2;Parent=t.m;exception=trans-splicing\n"
    )
    assert run_detect(doc, RepairContext(), [OP]) == []
