from ddbj_submission_tools.nuc.ddbj_gff import parse
from ddbj_submission_tools.nuc.ddbj_gff.repair import REGISTRY, RepairContext, run_detect, run_apply, apply_all_names

HDR = "##gff-version 3\n##sequence-region c 1 100000\n"
MARK = "rrna-partial-mark-location"
NOTE = "aligned only 58.21 percent of the 18S ribosomal RNA"


def _doc(body):
    return parse(HDR + body)


def _rrna(fid="r1", strand="+", extra="", ftype="rRNA", start=1000, end=2087):
    return (f"c\tpybarrnap:0.5.1\t{ftype}\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID={fid};product=18S ribosomal RNA (partial);Note={NOTE}{extra}\n")


def test_detect_finds_a_partial_marked_rrna():
    cands = run_detect(_doc(_rrna()), RepairContext(), [MARK])
    assert len(cands) == 1
    assert cands[0].feature_id == "r1"
    assert cands[0].payload["marker"] == NOTE
    assert "r1" in cands[0].detail


def test_detect_ignores_an_rrna_without_the_marker():
    doc = _doc("c\tpybarrnap:0.5.1\trRNA\t1000\t2600\t.\t+\t.\t"
               "ID=r2;product=18S ribosomal RNA\n")
    assert run_detect(doc, RepairContext(), [MARK]) == []


def test_detect_ignores_an_already_partial_rrna():
    doc = _doc(_rrna(extra=";partial=true;start_range=.,1000;end_range=2087,."))
    assert run_detect(doc, RepairContext(), [MARK]) == []


def test_detect_ignores_non_rrna_types_carrying_the_marker():
    doc = _doc(_rrna(fid="t1", ftype="tRNA") + _rrna(fid="m1", ftype="mRNA"))
    assert run_detect(doc, RepairContext(), [MARK]) == []


def test_detect_ignores_a_feature_without_an_id():
    doc = _doc(f"c\tpybarrnap:0.5.1\trRNA\t1000\t2087\t.\t+\t.\tproduct=18S;Note={NOTE}\n")
    assert run_detect(doc, RepairContext(), [MARK]) == []


def test_detect_matches_the_marker_in_any_attribute():
    # the validator scans every qualifier, not just Note
    doc = _doc("c\tpybarrnap:0.5.1\trRNA\t1000\t2087\t.\t+\t.\t"
               f"ID=r3;product=18S ribosomal RNA;inference=alignment:{NOTE}\n")
    assert [c.feature_id for c in run_detect(doc, RepairContext(), [MARK])] == ["r3"]


def test_apply_marks_both_ends_partial_plus_strand():
    doc = _doc(_rrna(strand="+"))
    changes = run_apply(doc, RepairContext(), [MARK])
    a = doc.feature_index["r1"].attributes
    assert a["partial"] == ["true"]
    assert a["start_range"] == [".,1000"]
    assert a["end_range"] == ["2087,."]
    assert [c.action for c in changes] == ["mark-rrna-partial"]


def test_apply_marks_both_ends_partial_minus_strand():
    doc = _doc(_rrna(strand="-"))
    run_apply(doc, RepairContext(), [MARK])
    a = doc.feature_index["r1"].attributes
    assert a["partial"] == ["true"]
    assert a["start_range"] == [".,1000"]
    assert a["end_range"] == ["2087,."]


def test_apply_does_not_change_the_feature_graph():
    doc = _doc(_rrna() + "c\tx\tgene\t5000\t6000\t.\t+\t.\tID=g\n")
    before = ([f.id for f in doc.features], sorted(doc.feature_index),
              [r.id for r in doc.roots],
              {f.id: sorted(c.id for c in f.children) for f in doc.features})
    run_apply(doc, RepairContext(), [MARK])
    after = ([f.id for f in doc.features], sorted(doc.feature_index),
             [r.id for r in doc.roots],
             {f.id: sorted(c.id for c in f.children) for f in doc.features})
    assert before == after


def test_apply_is_idempotent():
    doc = _doc(_rrna())
    run_apply(doc, RepairContext(), [MARK])
    assert run_detect(doc, RepairContext(), [MARK]) == []
    assert run_apply(doc, RepairContext(), [MARK]) == []


def test_mark_operation_is_destructive_and_excluded_from_apply_all():
    assert REGISTRY[MARK].destructive is True
    assert REGISTRY[MARK].requires_sequence is False
    assert MARK not in apply_all_names()


DROP = "rrna-partial-drop"


def test_drop_removes_a_parentless_rrna():
    # the heterosigma shape: pybarrnap rRNAs have no gene parent
    doc = _doc(_rrna() + "c\tx\tgene\t5000\t6000\t.\t+\t.\tID=g\n")
    changes = run_apply(doc, RepairContext(), [DROP])
    assert doc.feature_index.get("r1") is None
    assert [f.id for f in doc.features] == ["g"]          # the unrelated gene survives
    assert [r.id for r in doc.roots] == ["g"]
    assert [c.action for c in changes] == ["drop-rrna-partial"]


def test_drop_removes_the_rrna_children_too():
    doc = _doc(_rrna()
               + f"c\tpybarrnap:0.5.1\texon\t1000\t2087\t.\t+\t.\tID=r1.e1;Parent=r1\n")
    run_apply(doc, RepairContext(), [DROP])
    assert doc.feature_index.get("r1") is None
    assert doc.feature_index.get("r1.e1") is None
    assert [f.id for f in doc.features] == []


def test_drop_removes_a_gene_left_childless():
    doc = _doc("c\tx\tgene\t1000\t2087\t.\t+\t.\tID=rg\n"
               + _rrna(extra=";Parent=rg"))
    run_apply(doc, RepairContext(), [DROP])
    assert doc.feature_index.get("r1") is None
    assert doc.feature_index.get("rg") is None            # gene had only this rRNA
    assert [f.id for f in doc.features] == []


def test_drop_keeps_a_gene_with_a_surviving_sibling():
    doc = _doc("c\tx\tgene\t1000\t3000\t.\t+\t.\tID=rg\n"
               + _rrna(extra=";Parent=rg")
               + "c\tpybarrnap:0.5.1\trRNA\t2500\t3000\t.\t+\t.\t"
                 "ID=r_ok;product=5S ribosomal RNA;Parent=rg\n")
    run_apply(doc, RepairContext(), [DROP])
    assert doc.feature_index.get("r1") is None
    assert doc.feature_index.get("r_ok") is not None      # no marker -> untouched
    gene = doc.feature_index.get("rg")
    assert gene is not None and [c.id for c in gene.children] == ["r_ok"]


def test_drop_keeps_a_child_shared_with_a_surviving_parent():
    doc = _doc(_rrna(fid="r1")
               + "c\tpybarrnap:0.5.1\trRNA\t1000\t2087\t.\t+\t.\t"
                 "ID=r_ok;product=5S ribosomal RNA\n"
               + "c\tpybarrnap:0.5.1\texon\t1000\t2087\t.\t+\t.\tID=shared;Parent=r1,r_ok\n")
    run_apply(doc, RepairContext(), [DROP])
    assert doc.feature_index.get("r1") is None
    shared = doc.feature_index.get("shared")
    assert shared is not None and shared.parent_ids == ["r_ok"]
    assert shared.attributes["Parent"] == ["r_ok"]


def test_drop_is_idempotent():
    doc = _doc(_rrna())
    run_apply(doc, RepairContext(), [DROP])
    assert run_detect(doc, RepairContext(), [DROP]) == []
    assert run_apply(doc, RepairContext(), [DROP]) == []


def test_drop_operation_is_destructive_and_excluded_from_apply_all():
    assert REGISTRY[DROP].destructive is True
    assert REGISTRY[DROP].requires_sequence is False
    assert DROP not in apply_all_names()
