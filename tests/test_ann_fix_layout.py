from ddbj_submission_tools.nuc.ann_fix.annfile import AnnFile, _split
from ddbj_submission_tools.nuc.ann_fix.layout import collapse_repeats, merge_previous


def _ann(text: str) -> AnnFile:
    return AnnFile([_split(l) for l in text.splitlines(keepends=True)])


def _text(ann: AnnFile) -> str:
    return "".join(r.text() + r.eol for r in ann.rows)


def _t(*rows):
    return "".join("\t".join(r) + "\n" for r in rows)


def test_header_row_removed():
    ann = _ann(_t(("Entry", "Feature", "Location", "Qualifier", "Value"), ("COMMON", "DATE", "", "hold_date", "20270101")))
    collapse_repeats(ann)
    assert _text(ann) == _t(("COMMON", "DATE", "", "hold_date", "20270101"))


def test_keyword_written_twice_becomes_one_feature():
    ann = _ann(_t(("COMMON", "DATATYPE", "", "type", "WGS"), ("", "KEYWORD", "", "keyword", "WGS"),
                  ("", "KEYWORD", "", "keyword", "STANDARD_DRAFT"), ("", "DBLINK", "", "project", "PRJDB1")))
    collapse_repeats(ann)
    assert [f.key for f in ann.features()] == ["DATATYPE", "KEYWORD", "DBLINK"]
    assert ann.features()[1].qualifiers(ann.rows) == [("keyword", "WGS"), ("keyword", "STANDARD_DRAFT")]


def test_entry_and_feature_names_on_every_row_collapse():
    ann = _ann(_t(("COMMON", "SUBMITTER", "", "ab_name", "Mishima,T."), ("COMMON", "SUBMITTER", "", "ab_name", "Mishima,H."),
                  ("COMMON", "SUBMITTER", "", "contact", "Taro Mishima"), ("COMMON", "DATE", "", "hold_date", "20270101"),
                  ("seq1", "source", "1..100", "organism", "Escherichia coli"), ("seq1", "source", "", "mol_type", "genomic DNA")))
    collapse_repeats(ann)
    feats = ann.features()
    assert [(f.entry, f.key) for f in feats] == [("COMMON", "SUBMITTER"), ("COMMON", "DATE"), ("seq1", "source"), ("seq1", "source")]
    assert len(feats[0].qualifiers(ann.rows)) == 3
    assert [r.cols[0] for r in ann.rows] == ["COMMON", "", "", "", "seq1", ""]
    # a non-singleton feature (source) is left for a person: only the repeated entry name is cleared


def test_separated_duplicate_entry_is_reported_not_changed():
    text = _t(("COMMON", "DATE", "", "hold_date", "20270101"), ("seq1", "source", "1..100", "organism", "E. coli"),
              ("COMMON", "DBLINK", "", "project", "PRJDB1"))
    ann = _ann(text)
    log, review = collapse_repeats(ann)
    assert _text(ann) == text
    assert any("COMMON" in r and "line 3" in r for r in review), review


def test_split_reference_is_merged():
    ann = _ann(_t(("COMMON", "REFERENCE", "", "title", "A genome"), ("", "REFERENCE", "", "status", "Unpublished"),
                  ("", "REFERENCE", "", "ab_name", "Mishima,T."), ("", "REFERENCE", "", "ab_name", "Mishima,H."),
                  ("", "DATE", "", "hold_date", "20270101")))
    log, review = collapse_repeats(ann)
    assert [f.key for f in ann.features()] == ["REFERENCE", "DATE"] and review == []
    assert len(ann.features()[0].qualifiers(ann.rows)) == 4


def test_identical_references_keep_one():
    ref = [("", "REFERENCE", "", "title", "A genome"), ("", "", "", "status", "Unpublished"), ("", "", "", "ab_name", "Mishima,T.")]
    ann = _ann(_t(("COMMON", "DATE", "", "hold_date", "20270101"), *ref, *ref))
    log, review = collapse_repeats(ann)
    assert [f.key for f in ann.features()] == ["DATE", "REFERENCE"] and review == []


def test_two_complete_references_are_kept():
    text = _t(("COMMON", "REFERENCE", "", "title", "A genome"), ("", "", "", "status", "Unpublished"), ("", "", "", "ab_name", "Mishima,T."),
              ("", "REFERENCE", "", "title", "Another paper"), ("", "", "", "status", "Published"), ("", "", "", "ab_name", "Mishima,H."),
              ("", "", "", "year", "2025"))
    ann = _ann(text)
    log, review = collapse_repeats(ann)
    assert _text(ann) == text and review == []


def test_undecidable_references_are_reported():
    text = _t(("COMMON", "REFERENCE", "", "title", "A genome"), ("", "", "", "status", "Unpublished"), ("", "", "", "ab_name", "Mishima,T."),
              ("", "REFERENCE", "", "title", "A Genome."), ("", "", "", "ab_name", "Mishima,T."))
    ann = _ann(text)
    log, review = collapse_repeats(ann)
    assert _text(ann) == text
    assert any("REFERENCE" in r and "lines 1 and 4" in r for r in review), review


def test_merge_previous():
    ann = _ann(_t(("COMMON", "REFERENCE", "", "title", "A genome"), ("", "", "", "status", "Unpublished"),
                  ("", "REFERENCE", "", "ab_name", "Mishima,T.")))
    merge_previous(ann, 3)
    assert [f.key for f in ann.features()] == ["REFERENCE"]


def test_merge_previous_needs_the_same_key_right_before():
    import pytest
    ann = _ann(_t(("COMMON", "DATE", "", "hold_date", "20270101"), ("", "REFERENCE", "", "title", "A genome")))
    with pytest.raises(ValueError, match="REFERENCE"):
        merge_previous(ann, 2)
    with pytest.raises(ValueError, match="first row of a feature"):
        merge_previous(ann, 5)
