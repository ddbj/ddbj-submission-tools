"""nuc.organelle_ref.align: blastn best hits and position mapping (needs blastn on PATH)."""
import random

import pytest

from ddbj_submission_tools.nuc.organelle_ref import align

pytestmark = pytest.mark.skipif(not align.blast_available(), reason="blastn not on PATH")

_rng = random.Random(1)
REF = "".join(_rng.choice("ACGT") for _ in range(3000))


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def test_plus_strand_identity_and_coverage(tmp_path):
    hits = align.best_hits({"q": REF[1000:1600]}, {"R": REF}, tmp_path)
    h = hits["q"]["R"]
    assert (h.identity, h.coverage, h.strand) == (100.0, 1.0, 1)
    assert align.q2s(h, 1) == 1001 and align.s2q(h, 1600) == 600 and align.s2q(h, 50) is None


def test_minus_strand_maps_positions(tmp_path):
    h = align.best_hits({"q": revcomp(REF[1000:1600])}, {"R": REF}, tmp_path)["q"]["R"]
    assert h.strand == -1
    assert align.q2s(h, 1) == 1600 and align.q2s(h, 600) == 1001
    assert align.s2q(h, 1600) == 1 and align.s2q(h, 1001) == 600


def test_gap_in_query_maps_to_none(tmp_path):
    q = REF[1000:1300] + REF[1303:1600]            # 3-base deletion after ref position 1300
    h = align.best_hits({"q": q}, {"R": REF}, tmp_path)["q"]["R"]
    assert align.s2q(h, 1302) is None and align.s2q(h, 1304) == 301


def test_best_hit_per_reference_and_no_hit(tmp_path):
    rng2 = random.Random(2)
    other = "".join(rng2.choice("ACGT") for _ in range(3000))
    hits = align.best_hits({"q": REF[1000:1600], "z": other[:500]}, {"R": REF, "S": other}, tmp_path)
    assert set(hits["q"]) == {"R"} and set(hits["z"]) == {"S"}
