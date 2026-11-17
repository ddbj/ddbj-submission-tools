from ddbj_gff import parse
from ddbj_gff.repair import REGISTRY, RepairContext, run_detect, run_apply, apply_all_names

OP = "add-inference"
HDR = "##gff-version 3\n##sequence-region c 1 100000\n"

BODY = (
    "c\tAUGUSTUS\tgene\t1000\t2000\t.\t+\t.\tID=g1\n"
    "c\tAUGUSTUS\tmRNA\t1000\t2000\t.\t+\t.\tID=m1;Parent=g1\n"
    "c\tAUGUSTUS\tCDS\t1000\t2000\t.\t+\t0\tID=c1;Parent=m1\n"
    "c\ttRNAscan-SE\ttRNA\t5000\t5070\t.\t+\t.\tID=t1\n"
)


def _doc(body=BODY):
    return parse(HDR + body)


def _ids(cands):
    return sorted(c.feature_id for c in cands)


def test_default_filters_target_mrna_only():
    assert _ids(run_detect(_doc(), RepairContext(), [OP])) == ["m1"]


def test_apply_writes_the_default_value_on_the_mrna_only():
    doc = _doc()
    changes = run_apply(doc, RepairContext(), [OP])
    assert doc.feature_index["m1"].attributes["inference"] == ["ab initio prediction"]
    assert "inference" not in doc.feature_index["c1"].attributes
    assert "inference" not in doc.feature_index["t1"].attributes
    assert [c.action for c in changes] == ["add-inference"]


def test_inference_types_selects_other_types():
    ctx = RepairContext(inference_types="tRNA,ncRNA")
    assert _ids(run_detect(_doc(), ctx, [OP])) == ["t1"]


def test_inference_types_star_matches_every_id_bearing_feature():
    ctx = RepairContext(inference_types="*")
    assert _ids(run_detect(_doc(), ctx, [OP])) == ["c1", "g1", "m1", "t1"]


def test_a_star_entry_anywhere_in_the_list_matches_any_type():
    ctx = RepairContext(inference_types="mRNA,*")
    assert _ids(run_detect(_doc(), ctx, [OP])) == ["c1", "g1", "m1", "t1"]


def test_filter_lists_tolerate_whitespace_and_empty_entries():
    ctx = RepairContext(inference_types=" tRNA , , mRNA ")
    assert _ids(run_detect(_doc(), ctx, [OP])) == ["m1", "t1"]


def test_inference_sources_filters_by_gff_column_two():
    ctx = RepairContext(inference_types="*", inference_sources="AUGUSTUS")
    assert _ids(run_detect(_doc(), ctx, [OP])) == ["c1", "g1", "m1"]   # t1 is tRNAscan-SE


def test_an_empty_source_list_matches_any_source():
    ctx = RepairContext(inference_types="*", inference_sources="")
    assert _ids(run_detect(_doc(), ctx, [OP])) == ["c1", "g1", "m1", "t1"]


def test_an_existing_inference_is_never_overwritten():
    body = ("c\tAUGUSTUS\tmRNA\t1000\t2000\t.\t+\t.\t"
            "ID=m1;inference=similar to RNA sequence:INSD:LC123456.1\n")
    doc = _doc(body)
    assert run_detect(doc, RepairContext(), [OP]) == []
    assert run_apply(doc, RepairContext(), [OP]) == []
    assert doc.feature_index["m1"].attributes["inference"] == \
        ["similar to RNA sequence:INSD:LC123456.1"]


def test_a_custom_value_is_written_verbatim():
    ctx = RepairContext(inference_value="similar to AA sequence:INSD:BAA12345.1")
    doc = _doc()
    run_apply(doc, ctx, [OP])
    assert doc.feature_index["m1"].attributes["inference"] == \
        ["similar to AA sequence:INSD:BAA12345.1"]


def test_a_feature_without_an_id_is_skipped():
    doc = _doc("c\tAUGUSTUS\tmRNA\t1000\t2000\t.\t+\t.\tproduct=hypothetical protein\n")
    assert run_detect(doc, RepairContext(), [OP]) == []


def test_the_feature_graph_is_unchanged():
    doc = _doc()
    before = ([f.id for f in doc.features], sorted(doc.feature_index),
              [r.id for r in doc.roots],
              {f.id: sorted(c.id for c in f.children) for f in doc.features})
    run_apply(doc, RepairContext(inference_types="*"), [OP])
    after = ([f.id for f in doc.features], sorted(doc.feature_index),
             [r.id for r in doc.roots],
             {f.id: sorted(c.id for c in f.children) for f in doc.features})
    assert before == after


def test_apply_is_idempotent():
    doc = _doc()
    run_apply(doc, RepairContext(), [OP])
    assert run_detect(doc, RepairContext(), [OP]) == []
    assert run_apply(doc, RepairContext(), [OP]) == []


def test_operation_is_destructive_and_excluded_from_apply_all():
    assert REGISTRY[OP].destructive is True
    assert REGISTRY[OP].requires_sequence is False
    assert OP not in apply_all_names()


# --- detect_by_default: add-inference is opt-in, not part of a bare --detect ------

def test_bare_detect_returns_none_of_add_inference_but_keeps_a_default_op():
    # An mRNA with no UTR (exon and CDS share both ends) is a real defect for
    # utr-absent-to-partial-mrna (detect_by_default True) and, incidentally, also has no
    # /inference -- a bare sweep must report the former and omit the latter entirely.
    body = (
        "c\tAUGUSTUS\tgene\t1000\t2000\t.\t+\t.\tID=g1\n"
        "c\tAUGUSTUS\tmRNA\t1000\t2000\t.\t+\t.\tID=m1;Parent=g1\n"
        "c\tAUGUSTUS\texon\t1000\t2000\t.\t+\t.\tID=e1;Parent=m1\n"
        "c\tAUGUSTUS\tCDS\t1000\t2000\t.\t+\t0\tID=c1;Parent=m1\n"
    )
    doc = _doc(body)
    ops = {c.operation for c in run_detect(doc, RepairContext())}
    assert OP not in ops
    assert "utr-absent-to-partial-mrna" in ops


def test_registry_detect_by_default_flags():
    assert REGISTRY[OP].detect_by_default is False
    assert REGISTRY["utr-absent-to-partial-mrna"].detect_by_default is True
