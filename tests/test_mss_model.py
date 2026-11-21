from ddbj_submission_tools.nuc.gff2mss.model import MssQualifier, MssFeature


def test_model_construction():
    q = MssQualifier("locus_tag", "MPTK1_000010")
    f = MssFeature("CDS", "join(1..10,20..30)", [q])
    assert f.qualifiers[0].key == "locus_tag"
    assert f.location == "join(1..10,20..30)"


def test_defaults_independent():
    a = MssFeature("source", "1..9")
    b = MssFeature("source", "1..9")
    a.qualifiers.append(MssQualifier("organism", "X"))
    assert b.qualifiers == []
