"""Unit tests for _drop_subtree, the graph surgery shared by every drop operation.

The four operations that delete a feature (internal-utr-drop-mrna,
duplicate-mrna-keep-longest-cds, rrna-partial-drop, short-cds-drop-mrna) each own their
Change wording but share this helper. Their own test files cover the paths they exercise;
these tests pin the helper's contract directly, including two branches no operation test
reached: a feature with no gene parent, and remove_childless_gene=False.
"""

from ddbj_submission_tools.nuc.ddbj_gff import parse
from ddbj_submission_tools.nuc.ddbj_gff.repair.operations import _drop_subtree

HDR = "##gff-version 3\n##sequence-region c 1 100000\n"


def _doc(body):
    return parse(HDR + body)


GENE = "c\tX\tgene\t1000\t4000\t.\t+\t.\tID=g1\n"
M1 = "c\tX\tmRNA\t1000\t2000\t.\t+\t.\tID=m1;Parent=g1\n"
M2 = "c\tX\tmRNA\t3000\t4000\t.\t+\t.\tID=m2;Parent=g1\n"


def _ids(feats):
    return sorted(f.id for f in feats)


def test_a_parentless_feature_is_removed_with_no_gene_handling():
    # pybarrnap-style output: the RNA is a root feature, no gene layer at all
    doc = _doc("c\tX\trRNA\t1000\t2000\t.\t+\t.\tID=r1\n")
    exclusive, shared, dropped_gene = _drop_subtree(doc, doc.feature_index["r1"])
    assert (exclusive, shared, dropped_gene) == ([], [], None)
    assert doc.feature_index.get("r1") is None
    assert [f.id for f in doc.features] == []
    assert [r.id for r in doc.roots] == []


def test_descendants_with_no_outside_parent_are_all_deleted_and_reported():
    doc = _doc(GENE + M1
               + "c\tX\texon\t1000\t2000\t.\t+\t.\tID=m1.e1;Parent=m1\n"
               + "c\tX\tCDS\t1000\t2000\t.\t+\t0\tID=m1.c1;Parent=m1\n")
    exclusive, shared, dropped_gene = _drop_subtree(doc, doc.feature_index["m1"])
    assert _ids(exclusive) == ["m1.c1", "m1.e1"]
    assert shared == []
    assert dropped_gene is not None and dropped_gene.id == "g1"
    assert [f.id for f in doc.features] == []


def test_a_gene_with_a_surviving_sibling_is_kept_and_loses_only_the_dropped_child():
    doc = _doc(GENE + M1 + M2)
    exclusive, shared, dropped_gene = _drop_subtree(doc, doc.feature_index["m1"])
    assert dropped_gene is None
    gene = doc.feature_index.get("g1")
    assert gene is not None and [c.id for c in gene.children] == ["m2"]
    assert doc.feature_index.get("m1") is None


def test_remove_childless_gene_false_keeps_a_gene_that_would_otherwise_go():
    # what duplicate-mrna-keep-longest-cds relies on: the kept isoform is a sibling, and
    # even in a document where it is momentarily the only child the gene must survive
    doc = _doc(GENE + M1)
    exclusive, shared, dropped_gene = _drop_subtree(
        doc, doc.feature_index["m1"], remove_childless_gene=False)
    assert dropped_gene is None
    gene = doc.feature_index.get("g1")
    assert gene is not None and gene.children == []
    assert doc.feature_index.get("m1") is None


def test_a_child_shared_with_a_surviving_parent_survives_detached_on_all_three_links():
    doc = _doc(GENE + M1 + M2
               + "c\tX\texon\t1500\t1600\t.\t+\t.\tID=shared;Parent=m1,m2\n")
    exclusive, shared, dropped_gene = _drop_subtree(doc, doc.feature_index["m1"])
    assert exclusive == [] and _ids(shared) == ["shared"]
    ex = doc.feature_index.get("shared")
    assert ex is not None
    assert ex.parent_ids == ["m2"]                      # link 1
    assert ex.attributes["Parent"] == ["m2"]            # link 2
    assert [p.id for p in ex.parents] == ["m2"]         # link 3


def test_a_grandchild_below_a_shared_child_survives_with_it():
    # the shared exon survives, so the whole branch hanging off it must survive too --
    # the keep-set is grown to a fixed point rather than judged in one pass
    doc = _doc(GENE + M1 + M2
               + "c\tX\texon\t1500\t1600\t.\t+\t.\tID=shared;Parent=m1,m2\n"
               + "c\tX\tCDS\t1500\t1600\t.\t+\t0\tID=deep;Parent=shared\n")
    exclusive, shared, dropped_gene = _drop_subtree(doc, doc.feature_index["m1"])
    assert _ids(shared) == ["deep", "shared"]
    assert exclusive == []
    assert doc.feature_index.get("deep") is not None
    assert doc.feature_index["deep"].parent_ids == ["shared"]   # untouched, still valid


def test_a_whole_chain_below_a_shared_child_survives():
    # two levels below the rescued feature: a single extra relaxation pass would keep 'mid'
    # but still delete 'leaf'. Only a fixed point rescues the whole branch.
    doc = _doc(GENE + M1 + M2
               + "c\tX\texon\t1500\t1600\t.\t+\t.\tID=shared;Parent=m1,m2\n"
               + "c\tX\tCDS\t1500\t1600\t.\t+\t0\tID=mid;Parent=shared\n"
               + "c\tX\tCDS\t1520\t1560\t.\t+\t0\tID=leaf;Parent=mid\n")
    exclusive, shared, dropped_gene = _drop_subtree(doc, doc.feature_index["m1"])
    assert _ids(shared) == ["leaf", "mid", "shared"]
    assert exclusive == []
    assert all(doc.feature_index.get(i) is not None for i in ("shared", "mid", "leaf"))


def test_a_node_with_one_kept_and_one_deleted_parent_inside_the_subtree_survives():
    # both parents are inside the dropped subtree; 'shared' is rescued from outside and
    # 'doomed_e' is not, so the child must survive via 'shared' and lose only 'doomed_e'
    doc = _doc(GENE + M1 + M2
               + "c\tX\texon\t1500\t1600\t.\t+\t.\tID=shared;Parent=m1,m2\n"
               + "c\tX\texon\t1700\t1800\t.\t+\t.\tID=doomed_e;Parent=m1\n"
               + "c\tX\tCDS\t1500\t1800\t.\t+\t0\tID=kid;Parent=shared,doomed_e\n")
    exclusive, shared, dropped_gene = _drop_subtree(doc, doc.feature_index["m1"])
    assert _ids(exclusive) == ["doomed_e"] and _ids(shared) == ["kid", "shared"]
    kid = doc.feature_index.get("kid")
    assert kid is not None and kid.parent_ids == ["shared"]
    assert kid.attributes["Parent"] == ["shared"]


def test_a_shadow_features_id_does_not_strip_its_indexed_twins_link():
    """The parser tolerates one ID on two features when their types differ.

    It flags `id-type-mismatch` and keeps the second (shadow) out of `feature_index`, but
    still resolves parents through it. If the shadow lies in the dropped subtree, matching
    going parents by name alone would strip the *surviving* twin's name from a survivor's
    Parent, emitting an invalid empty `Parent=`.
    """
    doc = _doc(GENE + M1 + M2
               + "c\tX\texon\t3100\t3200\t.\t+\t.\tID=dup;Parent=m2\n"      # indexed, survives
               + "c\tX\tCDS\t1100\t1200\t.\t+\t0\tID=dup;Parent=m1\n"       # shadow, doomed
               + "c\tX\texon\t1500\t1600\t.\t+\t.\tID=e1;Parent=m1\n"
               + "c\tX\tCDS\t1500\t1600\t.\t+\t0\tID=kid;Parent=dup,e1\n")
    _drop_subtree(doc, doc.feature_index["m1"])
    kid = doc.feature_index.get("kid")
    assert kid is not None
    assert kid.parent_ids == ["dup"]                 # the surviving exon under m2
    assert kid.attributes["Parent"] == ["dup"]       # never an empty Parent=
    assert [p.id for p in kid.parents] == ["dup"]


def test_a_survivor_loses_its_link_to_a_deleted_intermediate_parent():
    # x survives via m2, but also named e1, which goes away with m1: leaving e1 in x's
    # Parent would write a reference to a feature no longer in the document
    doc = _doc(GENE + M1 + M2
               + "c\tX\texon\t1500\t1600\t.\t+\t.\tID=e1;Parent=m1\n"
               + "c\tX\tCDS\t1500\t1600\t.\t+\t0\tID=x;Parent=e1,m2\n")
    exclusive, shared, dropped_gene = _drop_subtree(doc, doc.feature_index["m1"])
    assert _ids(exclusive) == ["e1"] and _ids(shared) == ["x"]
    x = doc.feature_index.get("x")
    assert x is not None
    assert x.parent_ids == ["m2"]
    assert x.attributes["Parent"] == ["m2"]
    assert [p.id for p in x.parents] == ["m2"]
