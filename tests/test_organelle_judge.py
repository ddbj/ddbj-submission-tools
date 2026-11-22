"""nuc.organelle_ref.judge: identification verdicts, boundary comparison, fix.edit proposals."""
import random
from unittest.mock import ANY

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from ddbj_submission_tools.nuc.mssann.annotation import Feature, Qualifier
from ddbj_submission_tools.nuc.organelle_ref import align, judge
from ddbj_submission_tools.nuc.organelle_ref.align import Hsp
from ddbj_submission_tools.nuc.organelle_ref.refs import RefSummary

needs_blast = pytest.mark.skipif(not align.blast_available(), reason="blastn not on PATH")

_rng = random.Random(7)
REF = "".join(_rng.choice("ACGT") for _ in range(3000))


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def make_ref(acc="NC_000001", species="Apis cerana", extra=(), seq=REF, cds=True):
    rec = SeqRecord(Seq(seq), id=acc, name=acc, annotations={"organism": species, "molecule_type": "DNA"})
    rec.features += ([SeqFeature(FeatureLocation(1100, 1400, 1), type="CDS")] if cds else []) + [   # 1101..1400
        SeqFeature(FeatureLocation(500, 570, 1), type="tRNA"),                                       # 501..570
        SeqFeature(FeatureLocation(2000, 2300, 1), type="rRNA"), *extra]
    return RefSummary(acc, species, 1, ""), rec


def feat(key, loc, codon_start=None):
    q = [Qualifier("product", "x", 1)] + ([Qualifier("codon_start", str(codon_start), 1)] if codon_start else [])
    return Feature(key, loc, 1, q)


def hit(tmp_path, query, acc="NC_000001", ref=REF):
    return align.best_hits({"q": query}, {acc: ref}, tmp_path)["q"][acc]


def fake_hsp(identity, coverage=1.0, ref="R"):
    qlen = 600
    return Hsp("q", ref, identity, 500.0, 1, round(qlen * coverage), 1, round(qlen * coverage), qlen, "", "")


def results(findings):
    return [f.result for f in findings]


# -- identification ---------------------------------------------------------

def test_species_mismatch():
    hits = [(RefSummary("NC_2", "Apis b", 2, ""), fake_hsp(99.8)), (RefSummary("NC_1", "Apis a", 1, ""), fake_hsp(95.0))]
    assert "SPECIES_MISMATCH" in results(judge.identify("q", "Apis a", hits, judge.Thresholds()))


def test_other_species_match_when_claimed_absent():
    hits = [(RefSummary("NC_2", "Apis b", 2, ""), fake_hsp(99.5))]
    assert results(judge.identify("q", "Apis a", hits, judge.Thresholds())) == ["CLAIMED_ABSENT", "OTHER_SPECIES_MATCH"]


def test_claimed_below_coverage_is_not_called_absent():
    hits = [(RefSummary("NC_2", "Apis b", 2, ""), fake_hsp(99.3)), (RefSummary("NC_1", "Apis a", 1, ""), fake_hsp(100.0, 0.5))]
    assert "CLAIMED_ABSENT" not in results(judge.identify("q", "Apis a", hits, judge.Thresholds()))


def test_low_identity():
    hits = [(RefSummary("NC_1", "Apis a", 1, ""), fake_hsp(96.0)), (RefSummary("NC_2", "Apis b", 2, ""), fake_hsp(90.0))]
    assert results(judge.identify("q", "Apis a", hits, judge.Thresholds())) == ["LOW_IDENTITY"]


def test_claimed_best_gives_no_findings():
    hits = [(RefSummary("NC_1", "Apis a", 1, ""), fake_hsp(99.9)), (RefSummary("NC_2", "Apis b", 2, ""), fake_hsp(92.0))]
    assert judge.identify("q", "Apis a", hits, judge.Thresholds()) == []


def test_partial_map():
    hits = [(RefSummary("NC_1", "Apis a", 1, ""), fake_hsp(99.9, coverage=0.7))]
    assert "PARTIAL_MAP" in results(judge.identify("q", "Apis a", hits, judge.Thresholds()))


def test_no_hit_is_partial_map():
    assert results(judge.identify("q", "Apis a", [], judge.Thresholds())) == ["PARTIAL_MAP"]


def test_unresolved_species_gives_best_match_only():
    hits = [(RefSummary("NC_2", "Apis b", 2, ""), fake_hsp(100.0))]
    found = judge.identify("q", None, hits, judge.Thresholds())
    assert results(found) == ["BEST_MATCH"] and "Apis b" in found[0].detail
    assert "BEST_MATCH" not in judge.FINDING_RESULTS and "SPECIES_MISMATCH" in judge.FINDING_RESULTS


# -- boundaries -------------------------------------------------------------

Q = REF[1000:1600]                                          # ref 1001..1600 = query 1..600


@needs_blast
def test_matching_boundaries_give_no_findings(tmp_path):
    findings, rules = judge.boundary("q", [feat("CDS", "101..400")], Q, hit(tmp_path, Q), make_ref())
    assert (findings, rules) == ([], [])


@needs_blast
def test_shifted_start_gives_boundary_and_disabled_proposal(tmp_path):
    findings, rules = judge.boundary("q", [feat("CDS", "104..400")], Q, hit(tmp_path, Q), make_ref())
    assert results(findings) == ["BOUNDARY"]
    assert rules == [{"action": "replace_location", "name": "q CDS 104..400", "enabled": False, "entry": "q",
                      "feature": "CDS", "where_location": "^104\\.\\.400$", "pattern": "104..400",
                      "replacement": "101..400", "note": ANY}]


@needs_blast
def test_no_proposal_for_non_refseq_reference(tmp_path):
    h = hit(tmp_path, Q, acc="MN000001")
    findings, rules = judge.boundary("q", [feat("CDS", "104..400")], Q, h, make_ref(acc="MN000001"))
    assert results(findings) == ["BOUNDARY"] and rules == []


@needs_blast
def test_partial_end_ok_and_complete_end_inside(tmp_path):
    q = REF[1150:1600]                                      # ref 1151..1600; CDS 1101..1400 -> query <1..250
    findings, _ = judge.boundary("q", [feat("CDS", "<1..250", codon_start=2)], q, hit(tmp_path, q), make_ref())
    assert findings == []
    h = hit(tmp_path, Q)
    findings, _ = judge.boundary("q", [feat("tRNA", "<1..70")], Q, h, make_ref())     # tRNA 501..570 not here
    assert results(findings) == ["NOT_COMPARED"]
    findings, _ = judge.boundary("q", [feat("CDS", "<50..400")], Q, h, make_ref())
    assert "BOUNDARY" in results(findings)


def with_base(seq, pos, base):
    return seq[:pos - 1] + base + seq[pos:]


@needs_blast
def test_incomplete_stop_is_info_only_for_t_or_ta(tmp_path):
    # 101..398 is 298 nt: one base past the last full codon, at ref 1398
    for base, expected in (("T", "INCOMPLETE_STOP"), ("G", "BOUNDARY")):
        ref = with_base(REF, 1398, base)
        q = ref[1000:1600]
        findings, rules = judge.boundary("q", [feat("CDS", "101..398")], q, hit(tmp_path, q, ref=ref),
                                         make_ref(seq=ref))
        assert results(findings) == [expected], base
    assert "INCOMPLETE_STOP" not in judge.FINDING_RESULTS


@needs_blast
def test_cds_longer_than_reference_is_boundary(tmp_path):
    findings, _ = judge.boundary("q", [feat("CDS", "101..402")], Q, hit(tmp_path, Q), make_ref())
    assert results(findings) == ["BOUNDARY"]


@needs_blast
def test_codon_start_on_minus_strand(tmp_path):
    q = revcomp(REF[1150:1600])                             # query q = 1601 - ref position
    findings, rules = judge.boundary("q", [feat("CDS", "complement(201..>450)")], q, hit(tmp_path, q), make_ref())
    assert results(findings) == ["CODON_START"]
    assert rules == [{"action": "set_value", "name": "q CDS complement(201..>450)", "enabled": False, "entry": "q",
                      "feature": "CDS", "qualifier": "codon_start", "where_location": ANY, "value": "2",
                      "add_if_missing": True, "note": ANY}]


@needs_blast
def test_proposal_that_makes_5_end_partial_also_sets_codon_start(tmp_path):
    q = revcomp(REF[1150:1600])
    findings, rules = judge.boundary("q", [feat("CDS", "complement(201..447)")], q, hit(tmp_path, q), make_ref())
    assert "BOUNDARY" in results(findings) and "CODON_START" in results(findings)
    assert [(r["action"], r.get("replacement") or r.get("value")) for r in rules] == \
        [("set_value", "2"), ("replace_location", "complement(201..>450)")]      # set_value first: still matches


@needs_blast
def test_join_feature_compares_outer_ends_without_proposal(tmp_path):
    findings, rules = judge.boundary("q", [feat("CDS", "join(104..200,210..400)")], Q, hit(tmp_path, Q), make_ref())
    assert results(findings) == ["BOUNDARY"] and rules == []


@needs_blast
def test_reference_with_intron_is_not_compared_against_one_span(tmp_path):
    intron = SeqFeature(CompoundLocation([FeatureLocation(1100, 1200, 1), FeatureLocation(1700, 2000, 1)]), type="CDS")
    q = REF[1800:2400]                                      # ref 1801..2400; second exon 1701..2000 -> <1..200
    findings, rules = judge.boundary("q", [feat("CDS", "<1..200", codon_start=2)], q, hit(tmp_path, q),
                                     make_ref(cds=False, extra=[intron]))
    assert results(findings) == ["NOT_COMPARED"] and rules == []


@needs_blast
def test_trans_spliced_reference_is_not_compared(tmp_path):
    mixed = SeqFeature(CompoundLocation([FeatureLocation(200, 300, -1), FeatureLocation(1100, 1200, 1),
                                         FeatureLocation(2500, 2600, 1)]), type="CDS")
    findings, rules = judge.boundary("q", [feat("CDS", "101..200")], Q, hit(tmp_path, Q),
                                     make_ref(cds=False, extra=[mixed]))
    assert results(findings) == ["NOT_COMPARED"] and rules == []


@needs_blast
def test_feature_across_reference_origin_not_compared(tmp_path):
    across = SeqFeature(CompoundLocation([FeatureLocation(2950, 3000, 1), FeatureLocation(0, 50, 1)]), type="rRNA")
    q = REF[2900:3000]                                      # ref 2901..3000 = query 1..100
    findings, _ = judge.boundary("q", [feat("rRNA", "51..>100")], q, hit(tmp_path, q), make_ref(extra=[across]))
    assert results(findings) == ["NOT_COMPARED"] and "origin" in findings[0].detail


@needs_blast
def test_query_spanning_the_origin_maps_through_it(tmp_path):
    ext = REF + REF[:600]                                   # what the CLI aligns to for a circular genome
    q = REF[2700:] + REF[:300]                              # ref 2701..3000 + 1..300
    h = hit(tmp_path, q, ref=ext)
    assert h.coverage == 1.0
    findings, _ = judge.boundary("q", [feat("tRNA", "1..100")], q, h, make_ref())   # no tRNA there, but mapped
    assert results(findings) == ["NOT_COMPARED"] and "overlaps" in findings[0].detail
    findings, _ = judge.boundary("q", [feat("rRNA", "250..350")], q, h, make_ref())  # 2950..3000 + 1..50: wraps
    assert results(findings) == ["NOT_COMPARED"] and "origin" in findings[0].detail


@needs_blast
def test_reference_end_beyond_the_alignment_is_placed_by_reference_coordinates(tmp_path):
    rng = random.Random(99)
    q = "".join(rng.choice("ACGT") for _ in range(120)) + REF[1120:1600]   # query 121.. = ref 1121..; 1..120 diverged
    h = hit(tmp_path, q)
    assert h.qstart > 101                                                  # ref CDS start 1101 is not aligned
    findings, rules = judge.boundary("q", [feat("CDS", "125..400")], q, h, make_ref())
    assert results(findings) == ["BOUNDARY"] and "estimated" in findings[0].detail
    assert [r["replacement"] for r in rules] == ["101..400"]               # not <1..400


@needs_blast
def test_reference_end_in_a_query_gap_is_not_a_boundary(tmp_path):
    q = REF[1000:1100] + REF[1101:1600]                                    # ref 1101 (CDS start) deleted
    findings, rules = judge.boundary("q", [feat("CDS", "101..399")], q, hit(tmp_path, q), make_ref())
    assert (findings, rules) == ([], [])
