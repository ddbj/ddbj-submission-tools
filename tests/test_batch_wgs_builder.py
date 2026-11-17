"""batch_wgs_builder: output names, name collisions and the common JSON COMMENT."""

import json
import os

import pandas as pd
import pytest

from batch_wgs_builder.core import _output_prefix, build_sample_json, run

SEQ = "ACGT" * 50


def _row(values: dict[tuple[str, str], str]) -> pd.Series:
    return pd.Series(values)


def _write_inputs(tmp_path, rows: list[dict[str, str]], common: dict | None = None):
    columns = [("_", "_file_path"), ("DBLINK", "biosample"), ("source", "organism"), ("source", "strain")]
    lines = ["\t".join(c for c, _ in columns), "\t".join(q for _, q in columns)]
    for i, r in enumerate(rows):
        fasta = tmp_path / f"g{i}.fa"
        fasta.write_text(f">ctg1\n{SEQ}\n")
        lines.append("\t".join([str(fasta), r["biosample"], "Escherichia coli", r["strain"]]))
    tsv = tmp_path / "sample_list.tsv"
    tsv.write_text("\n".join(lines) + "\n")
    common_path = tmp_path / "common.json"
    common_path.write_text(json.dumps(common or {}))
    return str(tsv), str(common_path)


def test_output_prefix_keeps_only_safe_characters(tmp_path):
    cj = {"DBLINK": {"biosample": ["SAMD00000001"]}, "SOURCE": {"strain": "L2/434/Bu x:*"}}
    assert _output_prefix(str(tmp_path), cj) == os.path.join(str(tmp_path), "SAMD00000001_L2_434_Bu_x__")


def test_output_prefix_folds_non_ascii(tmp_path):
    cj = {"DBLINK": {"biosample": ["SAMD00000001"]}, "SOURCE": {"isolate": "株1"}}
    assert os.path.basename(_output_prefix(str(tmp_path), cj)) == "SAMD00000001__1"


@pytest.mark.parametrize("source_identifier,source,expected", [
    # SOURCE_IDENTIFIER に選んだ qualifier の値を使う
    ("cultivar", {"cultivar": "Koshihikari", "strain": "S1"}, "SAMD00000001_Koshihikari"),
    ("isolate", {"isolate": "I1", "strain": "S1"}, "SAMD00000001_I1"),
    # 選んだ qualifier の値が無い行は NO_IDENTIFIER（同じ名前の検査が効くように）
    ("cultivar", {"strain": "S1"}, "SAMD00000001_NO_IDENTIFIER"),
    # 未指定なら今までどおり strain、無ければ isolate
    ("", {"strain": "S1", "isolate": "I1"}, "SAMD00000001_S1"),
    (None, {"isolate": "I1"}, "SAMD00000001_I1"),
])
def test_output_prefix_uses_the_source_identifier(tmp_path, source_identifier, source, expected):
    cj = {"DBLINK": {"biosample": ["SAMD00000001"]}, "SOURCE": source}
    if source_identifier is not None:
        cj["SOURCE_IDENTIFIER"] = source_identifier
    assert os.path.basename(_output_prefix(str(tmp_path), cj)) == expected


def test_run_names_files_by_the_source_identifier(tmp_path):
    tsv, common = _write_inputs(tmp_path, [{"biosample": "SAMD00000001", "strain": "S1"}],
                                common={"SOURCE_IDENTIFIER": "strain"})
    lines = open(tsv).read().splitlines()
    lines[0] += "\tsource"
    lines[1] += "\tcultivar"
    lines[2] += "\tKoshihikari"
    open(tsv, "w").write("\n".join(lines) + "\n")
    common_json = json.loads(open(common).read())
    common_json["SOURCE_IDENTIFIER"] = "cultivar"
    open(common, "w").write(json.dumps(common_json))
    out = tmp_path / "out"
    run(tsv, common, str(out), submission_category="WGS")
    assert sorted(p.name for p in out.iterdir()) == ["SAMD00000001_Koshihikari.ann", "SAMD00000001_Koshihikari.fa"]


@pytest.mark.parametrize("common_comment,expected_head", [
    ([{"line": ["common 1", "x; y"]}], ["common 1", "x; y"]),
    ({"line": ["common 1"]}, ["common 1"]),
    ([{"line": "single line"}], ["single line"]),
])
def test_common_json_comment_comes_before_the_table_comment(common_comment, expected_head):
    row = _row({("_", "_file_path"): "g.fa", ("COMMENT", "line"): "a; b"})
    _, _, cj = build_sample_json(row, {"COMMENT": common_comment}, "WGS")
    assert cj["COMMENT"] == [{"line": expected_head + ["a", "b"]}]


def test_common_json_comment_without_table_comment():
    row = _row({("_", "_file_path"): "g.fa"})
    _, _, cj = build_sample_json(row, {"COMMENT": [{"line": ["only common"]}]}, "WGS")
    assert cj["COMMENT"] == [{"line": ["only common"]}]


def test_run_writes_safe_names(tmp_path):
    tsv, common = _write_inputs(tmp_path, [{"biosample": "SAMD00000001", "strain": "L2/434/Bu"}])
    out = tmp_path / "out"
    run(tsv, common, str(out), submission_category="WGS")
    assert sorted(os.listdir(out)) == ["SAMD00000001_L2_434_Bu.ann", "SAMD00000001_L2_434_Bu.fa"]


def test_run_stops_before_writing_when_names_collide(tmp_path):
    # "A/1" and "A:1" both become "A_1": the second would silently overwrite the first.
    tsv, common = _write_inputs(tmp_path, [
        {"biosample": "SAMD00000001", "strain": "A/1"},
        {"biosample": "SAMD00000001", "strain": "A:1"},
    ])
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="SAMD00000001_A_1"):
        run(tsv, common, str(out), submission_category="WGS")
    assert not out.exists() or os.listdir(out) == []
