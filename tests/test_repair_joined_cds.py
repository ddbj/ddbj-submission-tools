"""The completeness of a CDS exists across all its features, not within each one.

AGAT gives every CDS row its own ID, so one biological CDS is several Feature objects
sharing an mRNA parent. Judged apart, a fragment almost never begins with a start codon
and ends with a stop, so almost every fragment looks partial.
"""

from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.repair.context import RepairContext
from ddbj_gff.repair.registry import get_operation

HDR = "##gff-version 3\n##sequence-region s 1 100000\n"
OP = "missing-start-stop-to-partial-cds"

# Two CDS fragments, 1..6 and 11..16, with a 4 bp gap. Joined they read as one CDS;
# apart, fragment 1 has no stop and fragment 2 has no start.
#   plus-strand joined coding = seq[1..6] + seq[11..16]
COMPLETE_PLUS = Seq("ATGAAAGGGGGTTTAA")   # joined ATG AAA GTT TAA -> complete
NO_START_PLUS = Seq("CTTAAAGGGGGTTTAA")   # joined CTT AAA GTT TAA -> 5' partial only
NO_STOP_PLUS = Seq("ATGAAAGGGGGTTAAA")    # joined ATG AAA GTT AAA -> 3' partial only
#   minus-strand joined coding = revcomp(seq[11..16]) + revcomp(seq[1..6])
COMPLETE_MINUS = Seq("TTAAACGGGGTTTCAT")  # joined ATG AAA GTT TAA -> complete
NO_START_MINUS = Seq("TTAAACGGGGTTTAAG")  # joined CTT AAA GTT TAA -> 5' partial only


def _two_fragment(strand="+", extra1="", extra2=""):
    return HDR + "\n".join([
        f"s\tsrc\tgene\t1\t16\t.\t{strand}\t.\tID=g1",
        f"s\tsrc\tmRNA\t1\t16\t.\t{strand}\t.\tID=m1;Parent=g1",
        f"s\tsrc\tCDS\t1\t6\t.\t{strand}\t0\tID=cA;Parent=m1;transl_table=1{extra1}",
        f"s\tsrc\tCDS\t11\t16\t.\t{strand}\t0\tID=cB;Parent=m1;transl_table=1{extra2}",
    ]) + "\n"


def _detect(gff, seq):
    doc = parse(gff)
    ctx = RepairContext(sequences={"s": seq}, transl_table=1)
    return doc, ctx, get_operation(OP).detect(doc, ctx)


# --- the regression this work exists for ------------------------------------

def test_a_two_fragment_cds_complete_when_joined_is_not_partial():
    _, _, cands = _detect(_two_fragment(), COMPLETE_PLUS)
    assert cands == []


def test_the_same_on_the_minus_strand():
    _, _, cands = _detect(_two_fragment(strand="-"), COMPLETE_MINUS)
    assert cands == []


# --- which fragment carries the mark ----------------------------------------

def test_a_missing_start_marks_only_the_five_prime_fragment_plus_strand():
    doc, ctx, cands = _detect(_two_fragment(), NO_START_PLUS)
    assert [c.feature_id for c in cands] == ["cA"]        # genomically leftmost on +
    assert cands[0].payload["five"] is True
    assert cands[0].payload["three"] is False
    assert (cands[0].payload["start"], cands[0].payload["end"]) == (1, 6)   # cA's own extent
    get_operation(OP).apply(doc, ctx, None)
    assert doc.feature_index["cA"].attributes["start_range"] == [".,1"]
    assert "partial" not in doc.feature_index["cB"].attributes


def test_a_missing_stop_marks_only_the_three_prime_fragment_plus_strand():
    doc, ctx, cands = _detect(_two_fragment(), NO_STOP_PLUS)
    assert [c.feature_id for c in cands] == ["cB"]        # genomically rightmost on +
    assert cands[0].payload["three"] is True
    assert (cands[0].payload["start"], cands[0].payload["end"]) == (11, 16)
    get_operation(OP).apply(doc, ctx, None)
    assert doc.feature_index["cB"].attributes["end_range"] == ["16,."]
    assert "partial" not in doc.feature_index["cA"].attributes


def test_the_five_prime_fragment_is_the_genomically_rightmost_on_the_minus_strand():
    doc, ctx, cands = _detect(_two_fragment(strand="-"), NO_START_MINUS)
    assert [c.feature_id for c in cands] == ["cB"]        # 5' = genomic right on -
    assert cands[0].payload["five"] is True
    # 5' maps onto the genomic right end, so the mark is end_range
    get_operation(OP).apply(doc, ctx, None)
    assert doc.feature_index["cB"].attributes["end_range"] == ["16,."]
    assert "partial" not in doc.feature_index["cA"].attributes


# --- a chimeric group (disagreeing seqid or strand) is skipped, not guessed ------

def test_a_group_with_members_on_different_strands_yields_no_candidate():
    gff = HDR + "\n".join([
        "s\tsrc\tgene\t1\t16\t.\t+\t.\tID=g1",
        "s\tsrc\tmRNA\t1\t16\t.\t+\t.\tID=m1;Parent=g1",
        "s\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cA;Parent=m1;transl_table=1",
        "s\tsrc\tCDS\t11\t16\t.\t-\t0\tID=cB;Parent=m1;transl_table=1",
    ]) + "\n"
    _, _, cands = _detect(gff, NO_STOP_PLUS)
    assert cands == []


def test_a_group_with_members_on_different_seqids_yields_no_candidate():
    gff = (
        "##gff-version 3\n##sequence-region s 1 100000\n##sequence-region s2 1 100000\n"
        + "\n".join([
            "s\tsrc\tgene\t1\t16\t.\t+\t.\tID=g1",
            "s\tsrc\tmRNA\t1\t16\t.\t+\t.\tID=m1;Parent=g1",
            "s\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cA;Parent=m1;transl_table=1",
            "s2\tsrc\tCDS\t11\t16\t.\t+\t0\tID=cB;Parent=m1;transl_table=1",
        ]) + "\n"
    )
    doc = parse(gff)
    # both seqids are present in ctx.sequences -- the skip must come from the mismatch
    # itself, not from either seqid being absent.
    ctx = RepairContext(sequences={"s": NO_STOP_PLUS, "s2": NO_STOP_PLUS}, transl_table=1)
    cands = get_operation(OP).detect(doc, ctx)
    assert cands == []


# --- a three-fragment group: only the two ends are ever candidates ---------

# frag1 (1..6) = CTTAAA (no start), frag2 (11..16) = GGGGGG (filler, holds no end),
# frag3 (21..26) = GTTAAA (no stop). Joined = CTT AAA GGG GGG GTT AAA: 5' partial
# (no start), 3' partial (no stop); frag2's own bases never enter either check.
THREE_FRAGMENT_SEQ = Seq("CTTAAA" + "GGGG" + "GGGGGG" + "GGGG" + "GTTAAA")


def test_a_three_fragment_group_marks_only_the_two_ends():
    gff = HDR + "\n".join([
        "s\tsrc\tgene\t1\t26\t.\t+\t.\tID=g1",
        "s\tsrc\tmRNA\t1\t26\t.\t+\t.\tID=m1;Parent=g1",
        "s\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cA;Parent=m1;transl_table=1",
        "s\tsrc\tCDS\t11\t16\t.\t+\t0\tID=cMid;Parent=m1;transl_table=1",
        "s\tsrc\tCDS\t21\t26\t.\t+\t0\tID=cC;Parent=m1;transl_table=1",
    ]) + "\n"
    _, _, cands = _detect(gff, THREE_FRAGMENT_SEQ)
    by_id = {c.feature_id: c for c in cands}
    assert set(by_id) == {"cA", "cC"}                     # cMid untouched
    assert by_id["cA"].payload["five"] is True
    assert by_id["cA"].payload["three"] is False
    assert by_id["cC"].payload["five"] is False
    assert by_id["cC"].payload["three"] is True


# --- a group of one behaves exactly as before -------------------------------

def test_a_single_fragment_cds_still_gets_both_marks_on_itself():
    gff = HDR + "\n".join([
        "s\tsrc\tgene\t1\t12\t.\t+\t.\tID=g1",
        "s\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=m1;Parent=g1",
        "s\tsrc\tCDS\t1\t12\t.\t+\t0\tID=c1;Parent=m1;transl_table=1",
    ]) + "\n"
    doc, ctx, cands = _detect(gff, Seq("CTTAAAGTTAAA"))   # no start, no stop
    assert len(cands) == 1 and cands[0].feature_id == "c1"
    assert cands[0].payload["five"] is True and cands[0].payload["three"] is True
    get_operation(OP).apply(doc, ctx, None)
    a = doc.feature_index["c1"].attributes
    assert a["partial"] == ["true"] and a["start_range"] == [".,1"] and a["end_range"] == ["12,."]


def test_a_cds_with_no_mrna_parent_is_judged_alone():
    gff = HDR + "s\tsrc\tCDS\t1\t12\t.\t+\t0\tID=c1;transl_table=1\n"
    _, _, cands = _detect(gff, Seq("CTTAAAGTTAAA"))
    assert [c.feature_id for c in cands] == ["c1"]


# --- exemptions apply to the whole group ------------------------------------

def test_a_group_with_an_already_partial_member_is_skipped_whole():
    gff = _two_fragment(extra1=";partial=true;start_range=.%2C1")
    _, _, cands = _detect(gff, NO_STOP_PLUS)
    assert cands == []


def test_a_group_with_a_trans_spliced_member_is_skipped_whole():
    gff = _two_fragment(extra2=";exception=trans-splicing")
    _, _, cands = _detect(gff, NO_STOP_PLUS)
    assert cands == []


def test_codon_start_comes_from_the_five_prime_fragment_and_splits_the_marks():
    """phase on the 5'-most fragment governs the joined reading frame.

    Feature.codon_start reads ordered_spans()[0].phase, so a synthetic feature holding every
    span picks up the 5'-most fragment's phase automatically -- nothing needs to force it.
    phase 2 (codon_start 3) makes the first codon incomplete, and the frame shift also
    destroys the stop at the far end, so both ends are partial and each mark lands on its
    own terminal fragment.
    """
    gff = HDR + "\n".join([
        "s\tsrc\tgene\t1\t16\t.\t+\t.\tID=g1",
        "s\tsrc\tmRNA\t1\t16\t.\t+\t.\tID=m1;Parent=g1",
        "s\tsrc\tCDS\t1\t6\t.\t+\t2\tID=cA;Parent=m1;transl_table=1",
        "s\tsrc\tCDS\t11\t16\t.\t+\t0\tID=cB;Parent=m1;transl_table=1",
    ]) + "\n"
    doc, ctx, cands = _detect(gff, COMPLETE_PLUS)
    assert {c.feature_id: (c.payload["five"], c.payload["three"]) for c in cands} == {
        "cA": (True, False), "cB": (False, True)}
    get_operation(OP).apply(doc, ctx, None)
    assert doc.feature_index["cA"].attributes["start_range"] == [".,1"]
    assert "end_range" not in doc.feature_index["cA"].attributes
    assert doc.feature_index["cB"].attributes["end_range"] == ["16,."]
    assert "start_range" not in doc.feature_index["cB"].attributes


# --- a CDS shared by two mRNAs merges marks across groups, order-independently ---

# cA (1..6) and cShared (11..16) are mA's CDS children; cShared (11..16) and cB (21..26)
# are mB's. cShared is the 3'-terminal member of mA's group (rightmost of the pair) and
# the 5'-terminal member of mB's group (leftmost of the pair), so the two groups reach
# genuinely different, non-conflicting verdicts about it: mA says "no stop" (three), mB
# says "no start" (five). The merged result for cShared must be the union of both.
#   mA joined = cA + cShared = ATGAAA + GTTAAA -> ATG AAA GTT AAA: starts ATG, ends AAA (no stop)
#   mB joined = cShared + cB = GTTAAA + GTTTAA -> GTT AAA GTT TAA: starts GTT (no start), ends TAA
SHARED_SEQ = Seq("ATGAAAGGGGGTTAAAGGGGGTTTAA")  # cA=1-6 ATGAAA, cShared=11-16 GTTAAA, cB=21-26 GTTTAA


def _shared_cds_two_mrnas(swap=False):
    mrna_rows = [
        "s\tsrc\tmRNA\t1\t16\t.\t+\t.\tID=mA;Parent=g1",
        "s\tsrc\tmRNA\t11\t26\t.\t+\t.\tID=mB;Parent=g1",
    ]
    if swap:
        mrna_rows = mrna_rows[::-1]
    return HDR + "\n".join([
        "s\tsrc\tgene\t1\t26\t.\t+\t.\tID=g1",
        *mrna_rows,
        "s\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cA;Parent=mA;transl_table=1",
        "s\tsrc\tCDS\t11\t16\t.\t+\t0\tID=cShared;Parent=mA,mB;transl_table=1",
        "s\tsrc\tCDS\t21\t26\t.\t+\t0\tID=cB;Parent=mB;transl_table=1",
    ]) + "\n"


def test_a_cds_shared_by_two_mrnas_merges_marks_by_or_and_is_order_independent():
    def _summary(cands):
        return sorted((c.feature_id, c.payload["five"], c.payload["three"],
                       c.payload["start"], c.payload["end"]) for c in cands)

    _, _, cands_mA_first = _detect(_shared_cds_two_mrnas(swap=False), SHARED_SEQ)
    _, _, cands_mB_first = _detect(_shared_cds_two_mrnas(swap=True), SHARED_SEQ)
    expected = [("cShared", True, True, 11, 16)]   # union: five from mB, three from mA
    assert _summary(cands_mA_first) == expected
    assert _summary(cands_mB_first) == expected
