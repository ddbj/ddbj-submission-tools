from Bio.Seq import Seq
from ddbj_submission_tools.nuc.ddbj_gff import parse
from ddbj_submission_tools.nuc.ddbj_gff.repair import REGISTRY, RepairContext, run_detect, run_apply, apply_all_names

OP = "short-cds-drop-mrna"
HDR = "##gff-version 3\n##sequence-region c 1 100000\n"


def _doc(body):
    return parse(HDR + body)


def _tx(mid, cds_segments, gene="g1", extra_mrna="", extra_cds=""):
    """A gene-less mRNA row plus one CDS row per (start, end) segment."""
    lo = min(s for s, _ in cds_segments)
    hi = max(e for _, e in cds_segments)
    rows = (f"c\tAUGUSTUS\tmRNA\t{lo}\t{hi}\t.\t+\t.\t"
            f"ID={mid};Parent={gene}{extra_mrna}\n")
    for i, (s, e) in enumerate(cds_segments, 1):
        rows += (f"c\tAUGUSTUS\tCDS\t{s}\t{e}\t.\t+\t0\t"
                 f"ID={mid}.cds{i};Parent={mid}{extra_cds}\n")
    return rows


def _gene(gid="g1", lo=1000, hi=2000):
    return f"c\tAUGUSTUS\tgene\t{lo}\t{hi}\t.\t+\t.\tID={gid}\n"


def _ids(cands):
    return sorted(c.feature_id for c in cands)


# --- the threshold ----------------------------------------------------------

def test_a_cds_below_the_default_threshold_is_a_candidate():
    doc = _doc(_gene() + _tx("m1", [(1000, 1088)]))      # 89 bp
    cands = run_detect(doc, RepairContext(), [OP])
    assert _ids(cands) == ["m1"]
    assert cands[0].payload == {"cds_len": 89, "threshold": 90}


def test_a_cds_exactly_at_the_threshold_is_kept():
    doc = _doc(_gene() + _tx("m1", [(1000, 1089)]))      # 90 bp
    assert run_detect(doc, RepairContext(), [OP]) == []


def test_the_threshold_is_configurable():
    doc = _doc(_gene() + _tx("m1", [(1000, 1023)])                   # 24 bp -> dropped at 30
               + _gene("g2", 3000, 3032)
               + _tx("m2", [(3000, 3032)], gene="g2"))               # 33 bp -> kept at 30
    ctx = RepairContext(min_cds_length=30)
    assert _ids(run_detect(doc, ctx, [OP])) == ["m1"]


def test_multi_segment_cds_lengths_are_summed():
    # 30 + 30 = 60 bp total, below 90, though neither segment alone decides it
    doc = _doc(_gene() + _tx("m1", [(1000, 1029), (1500, 1529)]))
    cands = run_detect(doc, RepairContext(), [OP])
    assert cands[0].payload["cds_len"] == 60


# --- exemptions -------------------------------------------------------------

def test_an_mrna_whose_cds_is_already_partial_is_exempt():
    # ANN4300 applies to complete CDSs only
    doc = _doc(_gene() + _tx("m1", [(1000, 1023)], extra_cds=";partial=true;start_range=.%2C1000"))
    assert run_detect(doc, RepairContext(), [OP]) == []


def test_a_trans_spliced_mrna_is_exempt():
    # Feature.is_trans_spliced reads exception=trans-splicing, not a trans_splicing attribute
    doc = _doc(_gene() + _tx("m1", [(1000, 1023)], extra_mrna=";exception=trans-splicing"))
    assert run_detect(doc, RepairContext(), [OP]) == []


def test_an_mrna_with_a_trans_spliced_cds_child_is_exempt():
    doc = _doc(_gene() + _tx("m1", [(1000, 1023)], extra_cds=";exception=trans-splicing"))
    assert run_detect(doc, RepairContext(), [OP]) == []


def test_an_mrna_without_cds_children_is_not_a_candidate():
    doc = _doc(_gene() + "c\tAUGUSTUS\tmRNA\t1000\t2000\t.\t+\t.\tID=m1;Parent=g1\n")
    assert run_detect(doc, RepairContext(), [OP]) == []


def test_an_mrna_without_an_id_is_skipped():
    doc = _doc("c\tAUGUSTUS\tmRNA\t1000\t1023\t.\t+\t.\tproduct=hypothetical protein\n")
    assert run_detect(doc, RepairContext(), [OP]) == []


# --- apply ------------------------------------------------------------------

def test_apply_removes_the_transcript_its_children_and_a_childless_gene():
    doc = _doc(_gene() + _tx("m1", [(1000, 1023)])
               + "c\tAUGUSTUS\texon\t1000\t1023\t.\t+\t.\tID=m1.e1;Parent=m1\n")
    changes = run_apply(doc, RepairContext(), [OP])
    for gone in ("m1", "m1.cds1", "m1.e1", "g1"):
        assert doc.feature_index.get(gone) is None, gone
    assert [f.id for f in doc.features] == []
    assert [r.id for r in doc.roots] == []
    assert [c.action for c in changes] == ["drop-short-cds", "drop-short-cds"]


def test_apply_keeps_a_sibling_isoform_and_its_gene():
    doc = _doc(_gene(hi=4000) + _tx("m1", [(1000, 1023)])      # 24 bp -> dropped
               + _tx("m2", [(3000, 3200)]))                    # 201 bp -> kept
    run_apply(doc, RepairContext(), [OP])
    assert doc.feature_index.get("m1") is None
    assert doc.feature_index.get("m2") is not None
    gene = doc.feature_index.get("g1")
    assert gene is not None and [c.id for c in gene.children] == ["m2"]


def test_apply_keeps_a_child_shared_with_a_surviving_mrna():
    doc = _doc(_gene(hi=4000) + _tx("m1", [(1000, 1023)]) + _tx("m2", [(3000, 3200)])
               + "c\tAUGUSTUS\texon\t1000\t1023\t.\t+\t.\tID=shared;Parent=m1,m2\n")
    run_apply(doc, RepairContext(), [OP])
    shared = doc.feature_index.get("shared")
    assert shared is not None and shared.parent_ids == ["m2"]
    assert shared.attributes["Parent"] == ["m2"]


def test_apply_is_idempotent():
    doc = _doc(_gene() + _tx("m1", [(1000, 1023)]))
    run_apply(doc, RepairContext(), [OP])
    assert run_detect(doc, RepairContext(), [OP]) == []
    assert run_apply(doc, RepairContext(), [OP]) == []


def test_operation_is_destructive_and_excluded_from_apply_all():
    assert REGISTRY[OP].destructive is True
    assert REGISTRY[OP].requires_sequence is False
    assert OP not in apply_all_names()


# --- sequence-derived completeness (matches missing-start-stop-to-partial-cds) ----

def test_with_fasta_a_short_cds_with_start_and_stop_codons_is_still_a_candidate():
    # 24 bp, no partial= attribute: ATG ... TAA -- complete by the sequence test too.
    doc = _doc(_gene(lo=1, hi=24) + _tx("m1", [(1, 24)]))
    ctx = RepairContext(sequences={"c": Seq("ATGAAAGGGCCCTTTAAAGGGTAA")}, transl_table=1)
    assert _ids(run_detect(doc, ctx, [OP])) == ["m1"]


def test_with_fasta_a_short_cds_missing_a_stop_codon_is_not_a_candidate():
    # Same length and start codon, but the last codon is not a stop, and there is no
    # partial= attribute at all -- the sequence test must still exempt it.
    doc = _doc(_gene(lo=1, hi=24) + _tx("m1", [(1, 24)]))
    ctx = RepairContext(sequences={"c": Seq("ATGAAAGGGCCCTTTAAAGGGCCC")}, transl_table=1)
    assert run_detect(doc, ctx, [OP]) == []


def test_an_explicit_partial_attribute_wins_even_if_the_sequence_looks_complete():
    # partial=true stays an exemption even when start/stop codons are both present.
    doc = _doc(_gene(lo=1, hi=24)
               + _tx("m1", [(1, 24)], extra_cds=";partial=true;start_range=.%2C1"))
    ctx = RepairContext(sequences={"c": Seq("ATGAAAGGGCCCTTTAAAGGGTAA")}, transl_table=1)
    assert run_detect(doc, ctx, [OP]) == []


# --- fragmented CDS is judged as a whole (matches missing-start-stop-to-partial-cds) ----

# Two CDS fragments, 1..6 and 11..16, joined ATG AAA GTT TAA -- complete, and 12 bp total,
# far below the 90 bp threshold. Judged per fragment, fragment 1 has no stop and fragment 2
# no start, so both look partial and the transcript was wrongly exempted.
JOINED_SHORT = ("##gff-version 3\n##sequence-region s 1 100000\n"
                "s\tsrc\tgene\t1\t16\t.\t+\t.\tID=g1\n"
                "s\tsrc\tmRNA\t1\t16\t.\t+\t.\tID=m1;Parent=g1\n"
                "s\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cA;Parent=m1;transl_table=1\n"
                "s\tsrc\tCDS\t11\t16\t.\t+\t0\tID=cB;Parent=m1;transl_table=1\n")


def test_a_fragmented_but_complete_short_cds_is_a_candidate():
    doc = parse(JOINED_SHORT)
    ctx = RepairContext(sequences={"s": Seq("ATGAAAGGGGGTTTAA")}, transl_table=1)
    cands = run_detect(doc, ctx, [OP])
    assert [c.feature_id for c in cands] == ["m1"]
    assert cands[0].payload["cds_len"] == 12


def test_a_fragmented_incomplete_short_cds_is_still_exempt():
    # joined CTT AAA GTT TAA -- genuinely 5' partial, so ANN4300 would not fire on it
    doc = parse(JOINED_SHORT)
    ctx = RepairContext(sequences={"s": Seq("CTTAAAGGGGGTTTAA")}, transl_table=1)
    assert run_detect(doc, ctx, [OP]) == []


def test_an_explicit_partial_attribute_still_exempts_against_the_sequence():
    doc = parse(JOINED_SHORT.replace("ID=cA;Parent=m1", "ID=cA;Parent=m1;partial=true"))
    ctx = RepairContext(sequences={"s": Seq("ATGAAAGGGGGTTTAA")}, transl_table=1)
    assert run_detect(doc, ctx, [OP]) == []
