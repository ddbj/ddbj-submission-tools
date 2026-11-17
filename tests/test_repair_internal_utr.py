from ddbj_gff import parse
from ddbj_gff.model import cds_exon_intron_conflict
from ddbj_gff.repair import RepairContext, run_detect, run_apply, apply_all_names, DEFAULT_ORDER

HDR = "##gff-version 3\n##sequence-region c 1 100000\n"

# exons 50-200,300-400,500-650 ; CDS 100-200,500-600
# -> CDS intron (201,499) spans exon 300-400 => internal UTR
ANOM = (
    "c\tx\tgene\t50\t650\t.\t+\t.\tID=g\n"
    "c\tx\tmRNA\t50\t650\t.\t+\t.\tID=g.m;Parent=g\n"
    "c\tx\texon\t50\t200\t.\t+\t.\tID=g.e1;Parent=g.m\n"
    "c\tx\texon\t300\t400\t.\t+\t.\tID=g.e2;Parent=g.m\n"
    "c\tx\texon\t500\t650\t.\t+\t.\tID=g.e3;Parent=g.m\n"
    "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g.c1;Parent=g.m\n"
    "c\tx\tCDS\t500\t600\t.\t+\t0\tID=g.c2;Parent=g.m\n"
)

CLEAN = (
    "c\tx\tgene\t1000\t1600\t.\t+\t.\tID=h\n"
    "c\tx\tmRNA\t1000\t1600\t.\t+\t.\tID=h.m;Parent=h\n"
    "c\tx\texon\t1000\t1200\t.\t+\t.\tID=h.e1;Parent=h.m\n"
    "c\tx\texon\t1500\t1600\t.\t+\t.\tID=h.e2;Parent=h.m\n"
    "c\tx\tCDS\t1050\t1200\t.\t+\t0\tID=h.c1;Parent=h.m\n"
    "c\tx\tCDS\t1500\t1550\t.\t+\t0\tID=h.c2;Parent=h.m\n"
)

RECONCILE = "internal-utr-reconcile-exons"


def _doc(body):
    return parse(HDR + body)


def _exons(mrna):
    return sorted((s.start, s.end) for c in mrna.children if c.type == "exon" for s in c.spans)


def _cds(mrna):
    return sorted((s.start, s.end) for c in mrna.children if c.type == "CDS" for s in c.spans)


def test_detect_finds_anomalous_mrna():
    doc = _doc(ANOM + CLEAN)
    cands = run_detect(doc, RepairContext(), [RECONCILE])
    assert len(cands) == 1
    assert cands[0].feature_id == "g.m"
    assert cands[0].operation == RECONCILE
    assert cands[0].payload["conflicts"] == [[201, 499]]


def test_detect_silent_on_clean_mrna():
    doc = _doc(CLEAN)
    assert run_detect(doc, RepairContext(), [RECONCILE]) == []


def test_reconcile_rebuilds_exons_retaining_terminal_utr():
    doc = _doc(ANOM)
    changes = run_apply(doc, RepairContext(), [RECONCILE])
    m = doc.feature_index["g.m"]
    # interior exon 300-400 removed; terminal 5'UTR (50-99) and 3'UTR (601-650) kept
    assert _exons(m) == [(50, 200), (500, 650)]
    assert _cds(m) == [(100, 200), (500, 600)]          # CDS untouched
    assert cds_exon_intron_conflict(m) == []
    assert m.spans[0].start == 50 and m.spans[0].end == 650
    assert [c.action for c in changes] == ["reconcile-internal-utr"]


def test_reconcile_removes_old_exons_from_index():
    doc = _doc(ANOM)
    run_apply(doc, RepairContext(), [RECONCILE])
    for old in ("g.e1", "g.e2", "g.e3"):
        assert doc.feature_index.get(old) is None
    ids = [f.id for f in doc.features if f.type == "exon"]
    assert sorted(ids) == ["g.m.exon1", "g.m.exon2"]


def test_reconcile_is_idempotent():
    doc = _doc(ANOM)
    run_apply(doc, RepairContext(), [RECONCILE])
    assert run_detect(doc, RepairContext(), [RECONCILE]) == []


def test_reconcile_runs_before_internal_stop_in_apply_all():
    names = apply_all_names()
    assert names.index(RECONCILE) < names.index("internal-stop-to-misc")
    assert DEFAULT_ORDER.index(RECONCILE) < DEFAULT_ORDER.index("internal-stop-to-misc")


def test_reconcile_keeps_separate_terminal_utr_exons():
    # 5'UTR-only exon 10-50 and 3'UTR-only exon 850-900 are separate from the coding exons;
    # interior UTR exon 300-400 sits between coding exons.
    doc = _doc(
        "c\tx\tgene\t10\t900\t.\t+\t.\tID=k\n"
        "c\tx\tmRNA\t10\t900\t.\t+\t.\tID=k.m;Parent=k\n"
        "c\tx\texon\t10\t50\t.\t+\t.\tID=k.e0;Parent=k.m\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=k.e1;Parent=k.m\n"
        "c\tx\texon\t300\t400\t.\t+\t.\tID=k.e2;Parent=k.m\n"
        "c\tx\texon\t500\t600\t.\t+\t.\tID=k.e3;Parent=k.m\n"
        "c\tx\texon\t850\t900\t.\t+\t.\tID=k.e4;Parent=k.m\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=k.c1;Parent=k.m\n"
        "c\tx\tCDS\t500\t600\t.\t+\t0\tID=k.c2;Parent=k.m\n"
    )
    run_apply(doc, RepairContext(), [RECONCILE])
    m = doc.feature_index["k.m"]
    assert _exons(m) == [(10, 50), (100, 200), (500, 600), (850, 900)]
    assert cds_exon_intron_conflict(m) == []
    assert m.spans[0].start == 10 and m.spans[0].end == 900


def test_reconcile_exact_boundaries_no_terminal_utr():
    # every coding exon boundary equals the CDS boundary: no overhang to retain
    doc = _doc(
        "c\tx\tgene\t100\t600\t.\t+\t.\tID=n\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=n.m;Parent=n\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=n.e1;Parent=n.m\n"
        "c\tx\texon\t300\t400\t.\t+\t.\tID=n.e2;Parent=n.m\n"
        "c\tx\texon\t500\t600\t.\t+\t.\tID=n.e3;Parent=n.m\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=n.c1;Parent=n.m\n"
        "c\tx\tCDS\t500\t600\t.\t+\t0\tID=n.c2;Parent=n.m\n"
    )
    run_apply(doc, RepairContext(), [RECONCILE])
    m = doc.feature_index["n.m"]
    assert _exons(m) == _cds(m) == [(100, 200), (500, 600)]
    assert cds_exon_intron_conflict(m) == []


def test_reconcile_reports_dropped_introns():
    doc = _doc(ANOM
               + "c\tx\tintron\t201\t299\t.\t+\t.\tID=g.i1;Parent=g.m\n"
                 "c\tx\tintron\t401\t499\t.\t+\t.\tID=g.i2;Parent=g.m\n")
    changes = run_apply(doc, RepairContext(), [RECONCILE])
    assert len(changes) == 1
    assert "dropped 2 intron feature(s)" in changes[0].message
    assert doc.feature_index.get("g.i1") is None


DROP = "internal-utr-drop-mrna"

# a clean second isoform on the same gene g
CLEAN_ISOFORM = (
    "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.n;Parent=g\n"
    "c\tx\texon\t100\t200\t.\t+\t.\tID=g.n.e1;Parent=g.n\n"
    "c\tx\texon\t500\t600\t.\t+\t.\tID=g.n.e2;Parent=g.n\n"
    "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g.n.c1;Parent=g.n\n"
    "c\tx\tCDS\t500\t600\t.\t+\t0\tID=g.n.c2;Parent=g.n\n"
)


def test_drop_removes_mrna_children_and_childless_gene():
    doc = _doc(ANOM)
    changes = run_apply(doc, RepairContext(), [DROP])
    assert doc.feature_index.get("g.m") is None
    for child in ("g.e1", "g.e2", "g.e3", "g.c1", "g.c2"):
        assert doc.feature_index.get(child) is None
    assert doc.feature_index.get("g") is None            # gene left childless -> removed
    assert [f.id for f in doc.features] == []
    assert all(c.action == "drop-internal-utr" for c in changes)


def test_drop_keeps_clean_isoform_and_its_gene():
    doc = _doc(ANOM + CLEAN_ISOFORM)
    run_apply(doc, RepairContext(), [DROP])
    assert doc.feature_index.get("g.m") is None          # anomalous isoform dropped
    assert doc.feature_index.get("g.n") is not None      # clean isoform kept
    gene = doc.feature_index.get("g")
    assert gene is not None
    assert [c.id for c in gene.children] == ["g.n"]


def test_drop_is_idempotent():
    doc = _doc(ANOM)
    run_apply(doc, RepairContext(), [DROP])
    assert run_detect(doc, RepairContext(), [DROP]) == []


def test_drop_is_destructive_and_excluded_from_apply_all():
    from ddbj_gff.repair import REGISTRY
    assert REGISTRY[DROP].destructive is True
    assert DROP not in apply_all_names()


def test_apply_all_reconciles_rather_than_deleting():
    from ddbj_gff.repair import REGISTRY
    doc = _doc(ANOM)
    # run the non-sequence ops from the `all` list (the sequence ops would need a FASTA)
    names = [n for n in apply_all_names() if not REGISTRY[n].requires_sequence]
    run_apply(doc, RepairContext(), names)
    m = doc.feature_index.get("g.m")
    assert m is not None                                   # `all` never deletes annotation
    assert cds_exon_intron_conflict(m) == []                # it reconciled instead


def test_drop_keeps_child_shared_with_surviving_isoform():
    # g.e_shared is parented to BOTH the anomalous g.m and the clean g.n
    doc = _doc(
        "c\tx\tgene\t50\t650\t.\t+\t.\tID=g\n"
        "c\tx\tmRNA\t50\t650\t.\t+\t.\tID=g.m;Parent=g\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.n;Parent=g\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=g.e_shared;Parent=g.m,g.n\n"
        "c\tx\texon\t300\t400\t.\t+\t.\tID=g.e2;Parent=g.m\n"
        "c\tx\texon\t500\t650\t.\t+\t.\tID=g.e3;Parent=g.m\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g.c1;Parent=g.m\n"
        "c\tx\tCDS\t500\t600\t.\t+\t0\tID=g.c2;Parent=g.m\n"
        "c\tx\texon\t500\t600\t.\t+\t.\tID=g.n.e2;Parent=g.n\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g.n.c1;Parent=g.n\n"
        "c\tx\tCDS\t500\t600\t.\t+\t0\tID=g.n.c2;Parent=g.n\n"
    )
    changes = run_apply(doc, RepairContext(), [DROP])
    assert doc.feature_index.get("g.m") is None            # anomalous mRNA dropped
    assert doc.feature_index.get("g.e2") is None           # exclusive child dropped
    shared = doc.feature_index.get("g.e_shared")
    assert shared is not None                              # shared child survives
    assert shared.parent_ids == ["g.n"]                    # detached from the dropped mRNA
    assert shared.attributes["Parent"] == ["g.n"]
    assert all(id(p) != id(shared) or True for p in shared.parents)
    assert [p.id for p in shared.parents] == ["g.n"]
    assert any("kept 1 child feature(s) shared" in c.message for c in changes)


def test_internal_utr_fix_lives_in_repair_not_normalize():
    # the fix is a curation operation (repair), not always-safe canonicalization (normalize)
    # NB: import the submodule via importlib, not `from ddbj_gff.normalize import normalize`
    # -- the package's __init__.py does `from .normalize import normalize`, which rebinds
    # the package attribute to the *function*, shadowing the submodule of the same name.
    import importlib
    from ddbj_gff.normalize.config import NormalizeConfig
    norm_mod = importlib.import_module("ddbj_gff.normalize.normalize")
    assert not hasattr(NormalizeConfig(), "internal_utr_fix")
    assert "pass_fix_internal_utr" not in [p.__name__ for p in norm_mod.ALL_PASSES]
    assert not hasattr(norm_mod, "pass_fix_internal_utr")


def test_reconcile_keeps_child_shared_with_another_mrna():
    # g.e_shared belongs to both the anomalous g.m and the clean g.n
    doc = _doc(
        "c\tx\tgene\t50\t650\t.\t+\t.\tID=g\n"
        "c\tx\tmRNA\t50\t650\t.\t+\t.\tID=g.m;Parent=g\n"
        "c\tx\tmRNA\t100\t600\t.\t+\t.\tID=g.n;Parent=g\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=g.e_shared;Parent=g.m,g.n\n"
        "c\tx\texon\t300\t400\t.\t+\t.\tID=g.e2;Parent=g.m\n"
        "c\tx\texon\t500\t650\t.\t+\t.\tID=g.e3;Parent=g.m\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g.c1;Parent=g.m\n"
        "c\tx\tCDS\t500\t600\t.\t+\t0\tID=g.c2;Parent=g.m\n"
        "c\tx\texon\t500\t600\t.\t+\t.\tID=g.n.e2;Parent=g.n\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g.n.c1;Parent=g.n\n"
        "c\tx\tCDS\t500\t600\t.\t+\t0\tID=g.n.c2;Parent=g.n\n"
    )
    changes = run_apply(doc, RepairContext(), [RECONCILE])
    shared = doc.feature_index.get("g.e_shared")
    assert shared is not None                       # not deleted out from under g.n
    assert shared.parent_ids == ["g.n"]             # detached from the rebuilt mRNA only
    assert doc.feature_index.get("g.e2") is None    # exclusive child still removed
    assert any("kept 1 child feature(s) shared" in c.message for c in changes)
    # the clean isoform still has both its exons after the rebuild
    assert _exons(doc.feature_index["g.n"]) == [(100, 200), (500, 600)]


def test_drop_removes_grandchildren():
    doc = _doc(ANOM + "c\tx\trecoded_codon\t150\t152\t.\t+\t.\tID=g.rc;Parent=g.c1\n")
    run_apply(doc, RepairContext(), [DROP])
    assert doc.feature_index.get("g.c1") is None
    assert doc.feature_index.get("g.rc") is None     # grandchild goes with its parent
    assert [f.id for f in doc.features] == []


def test_trans_spliced_mrna_is_skipped():
    doc = _doc(
        "c\tx\tgene\t50\t650\t.\t+\t.\tID=t\n"
        "c\tx\tmRNA\t50\t650\t.\t+\t.\tID=t.m;Parent=t;exception=trans-splicing\n"
        "c\tx\texon\t50\t200\t.\t+\t.\tID=t.e1;Parent=t.m\n"
        "c\tx\texon\t300\t400\t.\t+\t.\tID=t.e2;Parent=t.m\n"
        "c\tx\texon\t500\t650\t.\t+\t.\tID=t.e3;Parent=t.m\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=t.c1;Parent=t.m\n"
        "c\tx\tCDS\t500\t600\t.\t+\t0\tID=t.c2;Parent=t.m\n"
    )
    assert run_detect(doc, RepairContext(), [RECONCILE]) == []
    assert run_detect(doc, RepairContext(), [DROP]) == []
