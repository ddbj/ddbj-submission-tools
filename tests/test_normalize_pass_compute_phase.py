from ddbj_submission_tools.nuc.ddbj_gff import parse
from ddbj_submission_tools.nuc.ddbj_gff.normalize.config import NormalizeConfig
from ddbj_submission_tools.nuc.ddbj_gff.normalize.normalize import normalize

HEAD = "##gff-version 3\nchr1\tL\tgene\t1\t999\t.\t{s}\t.\tID=g1\nchr1\tL\tmRNA\t1\t999\t.\t{s}\t.\tID=g1.t1;Parent=g1\n"


def _run(cds_rows, strand="+", on=True):
    gff = HEAD.format(s=strand) + "".join(
        f"chr1\tL\tCDS\t{a}\t{b}\t.\t{strand}\t{ph}\tID=g1.t1.cds{i};Parent=g1.t1\n"
        for i, (a, b, ph) in enumerate(cds_rows, 1))
    norm, rep = normalize(parse(gff), config=NormalizeConfig(taxid=1, compute_phase=on))
    phases = [norm.feature_index[f"g1.t1.cds{i}"].spans[0].phase for i in range(1, len(cds_rows) + 1)]
    return phases, [c for c in rep.applied if c.action == "set-phase"]


def test_plus_strand_phases():
    # lengths 54 (→ next 0), 65 (→ next (3-65%3)%3 = 1)
    phases, ch = _run([(97497, 97550, "."), (97600, 97664, "."), (97700, 97800, ".")])
    assert phases == [0, 0, 1] and len(ch) == 3


def test_minus_strand_walks_from_the_right():
    # transcript order on '-' is 300..310 (len 11 → next 1), 200..210, 100..110
    phases, ch = _run([(100, 110, "."), (200, 210, "."), (300, 310, ".")], strand="-")
    assert phases == [2, 1, 0]


def test_mixed_phase_is_untouched():
    phases, ch = _run([(1, 10, 2), (20, 30, ".")])
    assert phases == [2, None] and ch == []


def test_default_off():
    phases, ch = _run([(1, 10, "."), (20, 30, ".")], on=False)
    assert phases == [None, None] and ch == []
