from ddbj_gff import parse
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.normalize.normalize import normalize

HEAD = "##gff-version 3\nchr1\tA\tgene\t1\t30\t.\t{s}\t.\tID=g1\nchr1\tA\tmRNA\t1\t30\t.\t{s}\t.\tID=g1.t1;Parent=g1\n"


def _run_full(rows, strand="+", on=True):
    gff = HEAD.format(s=strand) + "".join(
        f"chr1\tA\t{t}\t{a}\t{b}\t.\t{strand}\t{ph}\tID={fid};Parent=g1.t1\n" for t, a, b, ph, fid in rows)
    norm, rep = normalize(parse(gff), config=NormalizeConfig(taxid=1, extend_cds_by_stop_codon=on))
    spans = sorted((s.start, s.end, s.phase) for s in norm.feature_index["g1.t1.cds"].spans)
    return spans, [c for c in rep.applied if c.action == "extend-cds"], rep


def _run(rows, strand="+", on=True):
    spans, ch, _ = _run_full(rows, strand, on)
    return spans, ch


def test_extends_adjacent_stop_on_plus():
    spans, ch = _run([("CDS", 1, 9, 0, "g1.t1.cds"), ("stop_codon", 10, 12, 0, "g1.t1.stop")])
    assert spans == [(1, 12, 0)] and len(ch) == 1 and ch[0].target == "g1.t1.cds"


def test_extends_adjacent_stop_on_minus():
    spans, ch = _run([("CDS", 10, 18, 0, "g1.t1.cds"), ("stop_codon", 7, 9, 0, "g1.t1.stop")], strand="-")
    assert spans == [(7, 18, 0)] and len(ch) == 1


def test_split_stop_codon_adds_span():
    spans, ch = _run([("CDS", 1, 9, 0, "g1.t1.cds"),
                      ("stop_codon", 10, 10, 0, "g1.t1.stop1"), ("stop_codon", 20, 21, 1, "g1.t1.stop2")])
    assert spans == [(1, 10, 0), (20, 21, 2)] and len(ch) == 1      # 10 nt precede the detached piece


def test_partially_overlapping_stop_is_not_appended():
    spans, ch, rep = _run_full([("CDS", 1, 11, 0, "g1.t1.cds"), ("stop_codon", 10, 12, 0, "g1.t1.stop")])
    assert spans == [(1, 11, 0)] and ch == []
    assert any(c.action == "needs-manual" and c.target == "g1.t1.cds" for c in rep.unresolved)


def test_single_detached_stop_is_not_appended():
    spans, ch, rep = _run_full([("CDS", 1, 9, 0, "g1.t1.cds"), ("stop_codon", 20, 22, 0, "g1.t1.stop")])
    assert spans == [(1, 9, 0)] and ch == []
    assert any(c.action == "needs-manual" for c in rep.unresolved)


def test_stop_inside_cds_is_noop():
    spans, ch = _run([("CDS", 1, 12, 0, "g1.t1.cds"), ("stop_codon", 10, 12, 0, "g1.t1.stop")])
    assert spans == [(1, 12, 0)] and ch == []


def test_default_off():
    spans, ch = _run([("CDS", 1, 9, 0, "g1.t1.cds"), ("stop_codon", 10, 12, 0, "g1.t1.stop")], on=False)
    assert spans == [(1, 9, 0)] and ch == []
