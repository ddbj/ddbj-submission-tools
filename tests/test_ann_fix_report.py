from pathlib import Path

import pytest

from ann_fix.report import Finding, group_by_code, parse_details

FX = Path(__file__).parent / "ann_fix_fixtures" / "details.txt"


def test_counts_per_code():
    groups = group_by_code(parse_details(FX.read_text()))
    assert {c: len(v) for c, v in groups.items()} == {
        "ANN5320": 2, "ANN4300": 2, "AXS6810": 2, "ANN2520": 1, "ANN2542": 1, "ANN0470": 1, "ANN2670": 1}
    assert list(groups)[:2] == ["ANN5320", "ANN4300"]        # first-appearance order


def test_fields_and_levels():
    f = group_by_code(parse_details(FX.read_text()))
    e = f["ANN5320"][0]
    assert e == Finding("ANN5320", "ERR", "NC_001142.9", "tRNA", 22138,
                        e.message) and e.message.startswith("The tRNA feature must not")
    assert f["ANN2520"][0].line == 1405 and f["ANN2520"][0].feature == "CDS"
    assert f["ANN2542"][0].line is None and f["ANN2542"][0].entry == "ALL"
    assert f["ANN0470"][0].line is None
    assert f["ANN2670"][0].level == "INF" and f["ANN2670"][0].line == 30


def test_headers_and_blank_lines_ignored():
    assert parse_details("=== Validation Details ===\nVersion: 1.0.70\n\n") == []


MULTI = """=== Validation Details (2 file sets) ===
Version: 1.0.95

1. alpha
[ ERROR ]
ANN5320:ERR:seq1:tRNA:Line[10]:The tRNA feature must not be completely contained within CDS exons. (Found: tRNA 5..80 in CDS 1..300)

2. beta
[ WARNING ]
ANN4300:WAR:seq9:CDS:Line[20]:Complete CDS must be at least 90 bases (30 aa) long unless supported by experiment or inference. (Length: 60 bases)
ANN0470:WAR:ALL:DBLINK:[FILE]:Missing SRA Run accession.
"""


def test_multi_file_report_records_file_set():
    fs = parse_details(MULTI)
    assert [(f.code, f.file_set) for f in fs] == [("ANN5320", "alpha"), ("ANN4300", "beta"), ("ANN0470", "beta")]


def test_multi_file_report_selects_one_file_set():
    assert [f.code for f in parse_details(MULTI, file_set="beta")] == ["ANN4300", "ANN0470"]


def test_unknown_file_set_names_the_sections():
    with pytest.raises(ValueError, match=r"gamma.*alpha, beta"):
        parse_details(MULTI, file_set="gamma")


def test_report_without_sections_keeps_all_for_any_file_set():
    # validator 1.0.70 writes no "N. NAME" header for a single file set
    assert len(parse_details(FX.read_text(), file_set="whatever")) == len(parse_details(FX.read_text()))
