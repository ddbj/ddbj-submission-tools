"""SUBMITTER: phone / fax / phext are not in the DDBJ annotation file spec
(https://www.ddbj.nig.ac.jp/ddbj/file-format.html#submitter), so create_common drops them."""
from ddbj_submission_tools.nuc.common.common_builder import create_common


def test_create_common_drops_phone_fax_phext_from_submitter():
    rows = create_common({"SUBMITTER": {"ab_name": ["Doe,J."], "contact": "J Doe", "phone": "+81-55",
                                        "fax": "1", "phext": "2", "email": "a@b.c"}})
    assert [r[3] for r in rows] == ["ab_name", "contact", "email"]
    assert rows[0][:2] == ["COMMON", "SUBMITTER"]


def test_create_common_submitter_with_only_phone_writes_nothing():
    assert create_common({"SUBMITTER": {"phone": "1"}}) == []


def test_create_common_keeps_phone_named_qualifiers_elsewhere():
    rows = create_common({"COMMENT": {"line": ["phone"]}, "DBLINK": {"project": "PRJDB1"}})
    assert [r[3] for r in rows] == ["line", "project"]
