import csv
from pathlib import Path

from ddbj_submission_tools.nuc.ann_fix.annfile import AnnFile
from ddbj_submission_tools.nuc.ann_meta.entries_tsv import write_entries
from ddbj_submission_tools.nuc.ann_meta.extract import extract, read_fasta_lengths
from ddbj_submission_tools.nuc.ann_meta.model import Entry, MetaDoc, Note

FIX = Path(__file__).parent / "ann_meta_fixtures"


def _doc(name: str, with_fasta: bool = True) -> MetaDoc:
    lengths = read_fasta_lengths(FIX / f"{name}.fasta") if with_fasta else None
    return extract(AnnFile.read(FIX / f"{name}.ann"), lengths)


def _ann_from_rows(tmp_path, rows: list[str]) -> AnnFile:
    p = tmp_path / "x.ann"
    p.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return AnnFile.read(p)


# -- fixtures -----------------------------------------------------------------

def test_dfast_draft():
    doc = _doc("dfast_draft")
    assert doc.category.value == "WGS"
    # organism / strain are empty placeholders, submitter_seqid is dropped,
    # ff_definition goes to each Entry
    assert doc.source == {"mol_type": "genomic DNA"}
    assert len(doc.entries) == 3
    assert [e.seq_id for e in doc.entries] == ["sequence1", "sequence2", "sequence3"]
    assert [e.type for e in doc.entries] == ["unplaced"] * 3
    assert [e.status for e in doc.entries] == ["partial"] * 3
    assert [e.length for e in doc.entries] == [3120, 2480, 1870]
    assert all(e.overrides == {} and e.extra == [] for e in doc.entries)
    assert doc.entries[0].ff_definition == "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"
    assert doc.source_identifier == "strain"
    assert doc.layout == "per_entry"
    assert doc.common["SUBMITTER"]["contact"] == ""
    assert doc.common["SUBMITTER"]["ab_name"] == ["", "", ""]
    assert doc.common["DBLINK"] == {"project": "", "sample": ""}
    assert doc.common["REFERENCE"] == [{"title": "", "ab_name": ["", "", ""], "status": "Unpublished"}]
    assert doc.common["COMMENT"] == [{"line": "Annotated by DFAST v.1.5.3 (https://dfast.ddbj.nig.ac.jp)"}]
    assert doc.common["ST_COMMENT"] == {"tagset_id": "Genome-Assembly-Data", "Assembly Method": "",
                                        "Genome Coverage": "", "Sequencing Technology": ""}
    assert doc.untouched_common_features == ["DATATYPE", "KEYWORD"]
    assert doc.common_extra == {}
    assert doc.unknown_qualifiers == []
    assert doc.notes == []


def test_dfast_complete():
    doc = _doc("dfast_complete")
    assert doc.category.value == "GNM"
    assert [e.type for e in doc.entries] == ["chromosome", "plasmid"]
    assert [e.seq_name for e in doc.entries] == ["1", "pXYZ"]
    assert [e.status for e in doc.entries] == ["complete", "complete"]
    assert [e.topology for e in doc.entries] == ["circular", "circular"]
    assert doc.source_identifier == "strain"
    assert doc.source == {"organism": "Escherichia coli", "strain": "NIG-0001", "mol_type": "genomic DNA",
                          "geo_loc_name": "Japan:Shizuoka, Mishima", "collection_date": "2025-06-12",
                          "host": "Homo sapiens", "isolation_source": "feces"}
    assert all(e.overrides == {} for e in doc.entries)
    assert doc.common["DBLINK"] == {"project": "PRJDB99901", "sample": "SAMD99000001",
                                    "sequence_read_archive": ["DRR990001"]}
    assert doc.common["REFERENCE"][0]["year"] == "2026"
    assert doc.untouched_common_features == []
    assert doc.notes == []


def test_marker_set():
    doc = _doc("marker_set")
    # marker_set.ann intentionally has no "organelle mitochondrion" on the COI
    # entries: with /organelle on every entry detect_category would say MISC.
    assert doc.category.value is None
    assert doc.source == {"mol_type": "genomic DNA"}
    for e in doc.entries:
        assert "organism" in e.overrides and "isolate" in e.overrides
        assert e.ff_definition == "@@[organism]@@ @@[isolate]@@ COI gene, partial cds"
        assert e.type == "unplaced" and e.status == "partial"
    e1, e2, e3 = doc.entries
    assert e1.overrides["organism"] == "Stenopsyche marmorata"
    assert e2.overrides["isolate"] == "SZ-014"
    assert e3.overrides["geo_loc_name"] == "Japan:Shizuoka, Kise River"
    # COI_02 has two note rows: kept in order in extra, not in overrides
    assert Note("repeated_qualifier", "COI_02", "note") in doc.notes
    assert e2.extra == [("note", "larva"), ("note", "voucher deposited in NSMT")]
    assert "note" not in e2.overrides
    assert e1.overrides["note"] == "voucher deposited in NSMT"
    assert e3.overrides["note"] == "adult male"
    assert e1.extra == [] and e3.extra == []
    assert doc.source_identifier == "isolate"
    assert doc.layout == "per_entry"
    assert "DBLINK" not in doc.common


def test_handmade_common_source_and_unknowns():
    # handmade.ann (and .fasta) use CRLF line endings.
    assert b"\r\n" in (FIX / "handmade.ann").read_bytes()
    doc = _doc("handmade")
    assert doc.layout == "common_source"
    assert doc.category.value == "WGS"
    assert "KEYWORD" in doc.untouched_common_features
    assert doc.common_extra["SUBMITTER"] == [("fax", "81-55-981-0000")]
    assert "DBLINK" not in doc.common
    assert doc.common["DATE"] == {"hold_date": "20271231"}
    assert doc.common["SUBMITTER"]["zip"] == "411-0000"  # no stray "\r"
    assert doc.source == {"organism": "Bacillus subtilis", "mol_type": "genomic DNA",
                          "strain": "HS-7", "sub_strain": "HS-7a"}
    assert doc.unknown_qualifiers == ["sub_strain"]
    assert [e.seq_id for e in doc.entries] == ["contig01", "contig02", "contig03"]
    assert [e.length for e in doc.entries] == [960, 720, 600]
    assert all(e.ff_definition == "@@[organism]@@ @@[strain]@@ DNA, contig: @@[entry]@@" for e in doc.entries)
    assert doc.source_identifier == "strain"
    assert doc.notes == []


def test_lengths_none_without_fasta():
    doc = _doc("dfast_draft", with_fasta=False)
    assert [e.length for e in doc.entries] == [None, None, None]


def test_read_fasta_lengths_crlf():
    assert read_fasta_lengths(FIX / "handmade.fasta") == {"contig01": 960, "contig02": 720, "contig03": 600}


# -- edge cases built in the test -------------------------------------------

def test_no_source_note(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "COMMON\tSUBMITTER\t\tab_name\tSuzuki,H.",
        "seq1\tCDS\t1..300\tproduct\thypothetical protein",
        "\t\t\ttransl_table\t11",
    ])
    doc = extract(ann)
    assert doc.notes[0].code == "no_source"
    assert doc.source_layout == "none"
    # the entries are still listed (with no source rows) so a caller can send them back
    assert [(e.seq_id, e.source_rows) for e in doc.entries] == [("seq1", [])]
    assert doc.source == {}


def test_no_source_lists_ann_and_fasta_entries(tmp_path):
    """No source anywhere: the entries of the ann and of the FASTA (FASTA order), with TOPOLOGY and DBLINK."""
    ann = _ann_from_rows(tmp_path, [
        "COMMON\tSUBMITTER\t\tab_name\tSuzuki,H.",
        "seq2\tTOPOLOGY\t\tcircular\t",
        "\tDBLINK\t\tbiosample\tSAMD1",
        "\tCDS\t1..300\tproduct\thypothetical protein",
        "seq1\tCDS\t1..90\tproduct\thypothetical protein",
    ])
    doc = extract(ann, {"seq1": 90, "seq2": 300, "seq3": 120})
    assert doc.notes[0].code == "no_source"
    assert doc.source_layout == "none"
    assert [e.seq_id for e in doc.entries] == ["seq1", "seq2", "seq3"]
    assert [e.source_rows for e in doc.entries] == [[], [], []]
    assert [e.topology for e in doc.entries] == ["", "circular", ""]
    assert [e.dblink for e in doc.entries] == [{}, {"sample": "SAMD1"}, {}]
    assert [e.length for e in doc.entries] == [90, 300, 120]


def test_both_sources_note(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "COMMON\tsource\t1..E\torganism\tBacillus subtilis",
        "\t\t\tmol_type\tgenomic DNA",
        "\t\t\tstrain\tHS-7",
        "seq1\tsource\t1..300\tstrain\tHS-8",
        "\tCDS\t1..300\tproduct\thypothetical protein",
        "seq2\tCDS\t1..300\tproduct\thypothetical protein",
    ])
    doc = extract(ann)
    assert Note("both_common_and_entry_source", None, "") in doc.notes
    assert doc.layout == "per_entry"
    assert doc.source == {"organism": "Bacillus subtilis", "mol_type": "genomic DNA"}
    assert [e.seq_id for e in doc.entries] == ["seq1", "seq2"]
    assert doc.entries[0].overrides == {"strain": "HS-8"}
    assert doc.entries[1].overrides == {"strain": "HS-7"}


def test_flag_qualifier_is_true_and_role_restoration(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "mt\tsource\t1..16000\torganism\tHomo sapiens",
        "\t\t\tmol_type\tgenomic DNA",
        "\t\t\torganelle\tmitochondrion",
        "\t\t\tenvironmental_sample\t",
        "\t\t\tsubmitter_seqid\tmt1",
        "\t\t\tff_definition\t@@[organism]@@ @@[cultivar]@@ mitochondrial DNA, complete genome",
    ])
    doc = extract(ann)
    e = doc.entries[0]
    assert (e.type, e.seq_name, e.status, e.topology) == ("organelle", "mitochondrion", "complete", "")
    assert doc.source == {"organism": "Homo sapiens", "mol_type": "genomic DNA", "environmental_sample": True}
    assert doc.source_identifier == "cultivar"
    assert e.length is None


def test_to_json_from_json_roundtrip():
    for name in ("dfast_draft", "dfast_complete", "marker_set", "handmade"):
        doc = _doc(name)
        d = doc.to_json()
        assert MetaDoc.from_json(d) == doc
        import json
        assert MetaDoc.from_json(json.loads(json.dumps(d))) == doc
    d = _doc("handmade").to_json()
    assert d["common_extra"] == {"SUBMITTER": [["fax", "81-55-981-0000"]]}
    assert d["unknown_qualifiers"] == ["sub_strain"]


def test_write_entries_tsv(tmp_path):
    doc = _doc("marker_set")
    out = tmp_path / "entries.tsv"
    write_entries(doc.entries, out)
    with open(out, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = rows[0]
    assert header[:6] == ["#seq_id", "type", "seq_name", "status", "topology", "ff_definition"]
    assert header[6:] == ["collection_date", "geo_loc_name", "isolate", "note", "note", "organism"]
    by_id = {r[0]: dict(zip(header, r)) for r in rows[1:]}
    coi1 = next(r for r in rows if r[0] == "COI_01")
    assert coi1[9:11] == ["voucher deposited in NSMT", ""]
    coi2 = next(r for r in rows if r[0] == "COI_02")
    assert coi2[9:11] == ["larva", "voucher deposited in NSMT"]    # one column per occurrence
    assert by_id["COI_02"]["organism"] == "Hydropsyche orientalis"
    assert by_id["COI_03"]["type"] == "unplaced"
    assert by_id["COI_03"]["ff_definition"] == "@@[organism]@@ @@[isolate]@@ COI gene, partial cds"


def test_write_entries_flag_and_none(tmp_path):
    out = tmp_path / "e.tsv"
    write_entries([Entry("s1", None, overrides={"environmental_sample": True}),
                   Entry("s2", 10, type="plasmid", seq_name="p1", status="complete", topology="circular",
                         extra=[("sub_strain", "a")])], out)
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "#seq_id\ttype\tseq_name\tstatus\ttopology\tff_definition\tenvironmental_sample\tsub_strain"
    assert lines[1] == "s1\tunplaced\t\tpartial\t\t\tyes\t"
    assert lines[2] == "s2\tplasmid\tp1\tcomplete\tcircular\t\t\ta"


# -- fix round 1: every dropped row gets a note ------------------------------

def test_unnamed_qualifier_row_is_noted(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "COMMON\tSUBMITTER\t\tab_name\tSuzuki,H.",
        "\t\t\t\tstray common value",
        "seq1\tsource\t1..300\torganism\tBacillus subtilis",
        "\t\t\t\torphan value",
        "\t\t\tmol_type\tgenomic DNA",
    ])
    doc = extract(ann)
    assert Note("unnamed_qualifier", None, "stray common value") in doc.notes
    assert Note("unnamed_qualifier", "seq1", "orphan value") in doc.notes
    assert doc.source == {"organism": "Bacillus subtilis", "mol_type": "genomic DNA"}


def test_topology_not_circular_is_noted(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "seq1\tsource\t1..300\torganism\tBacillus subtilis",
        "\tTOPOLOGY\t\tlinear",
    ])
    doc = extract(ann)
    assert doc.entries[0].topology == ""
    assert Note("topology_not_circular", "seq1", "linear") in doc.notes


def test_common_source_without_entries_is_noted(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "COMMON\tsource\t1..E\torganism\tBacillus subtilis",
        "\t\t\tmol_type\tgenomic DNA",
    ])
    doc = extract(ann)
    assert doc.entries == []
    assert Note("common_source_without_entries", None, "") in doc.notes


def test_repeated_common_feature_is_noted(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "COMMON\tSUBMITTER\t\tab_name\tSuzuki,H.",
        "\tSUBMITTER\t\tab_name\tTanaka,K.",
        "seq1\tsource\t1..300\torganism\tBacillus subtilis",
    ])
    doc = extract(ann)
    assert Note("repeated_common_feature", None, "SUBMITTER") in doc.notes
    assert doc.common["SUBMITTER"] == {"ab_name": ["Suzuki,H."]}
    assert doc.common_extra["SUBMITTER"] == [("ab_name", "Tanaka,K.")]


def test_multiple_roles_is_noted(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "seq1\tsource\t1..300\torganism\tBacillus subtilis",
        "\t\t\tchromosome\t1",
        "\t\t\tplasmid\tp1",
    ])
    doc = extract(ann)
    e = doc.entries[0]
    assert (e.type, e.seq_name) == ("chromosome", "1")
    assert e.extra == [("plasmid", "p1")]
    assert Note("multiple_roles", "seq1", "plasmid") in doc.notes


def test_entry_without_source_and_not_in_fasta(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "seq1\tsource\t1..300\torganism\tBacillus subtilis",
        "seq2\tCDS\t1..300\tproduct\thypothetical protein",
    ])
    doc = extract(ann, {"seq1": 300, "seq3": 120})
    assert Note("entry_without_source", "seq2", "") in doc.notes
    assert Note("entry_without_source", "seq3", "") in doc.notes
    assert Note("entry_not_in_fasta", "seq2", "") in doc.notes
    assert [e.seq_id for e in doc.entries] == ["seq1", "seq2", "seq3"]
    assert doc.source == {"organism": "Bacillus subtilis"}


def test_source_location_mismatch_is_noted(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "seq1\tsource\t1..290\torganism\tBacillus subtilis",
    ])
    doc = extract(ann, {"seq1": 300})
    assert Note("source_location", "seq1", "1..290") in doc.notes


def test_unknown_empty_valued_extra_survives_tsv(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "seq1\tsource\t1..300\torganism\tBacillus subtilis",
        "\t\t\tmy_flag\t",
        "seq2\tsource\t1..300\torganism\tBacillus subtilis",
    ])
    doc = extract(ann)
    assert doc.entries[0].extra == [("my_flag", True)]
    assert MetaDoc.from_json(doc.to_json()) == doc
    out = tmp_path / "e.tsv"
    write_entries(doc.entries, out)
    header, row1, row2 = out.read_text(encoding="utf-8").splitlines()
    assert header.split("\t")[-1] == "my_flag"
    assert row1.split("\t")[-1] == "yes"
    assert row2.split("\t")[-1] == ""


# -- final review fixes ------------------------------------------------------------

COMPLETE_NO_CHR = [
    "COMMON\tSUBMITTER\t\tab_name\tSuzuki,H.",
    "c1\tTOPOLOGY\t\tcircular\t",
    "\tsource\t1..5000\torganism\tBacillus subtilis",
    "\t\t\tmol_type\tgenomic DNA",
    "\t\t\tff_definition\t@@[organism]@@ DNA, complete genome",
    "\tCDS\t1..300\tproduct\thypothetical protein",
]


def test_unnamed_single_chromosome(tmp_path):
    # I3a: a complete genome without /chromosome is a chromosome with an empty name
    doc = extract(_ann_from_rows(tmp_path, COMPLETE_NO_CHR), {"c1": 5000})
    e = doc.entries[0]
    assert (e.type, e.seq_name, e.status, e.topology) == ("chromosome", "", "complete", "circular")


def test_unnamed_chromosome_by_ff_definition_or_topology(tmp_path):
    rows = [
        "c1\tsource\t1..500\torganism\tBacillus subtilis",
        "\t\t\tff_definition\t@@[organism]@@ DNA, chromosome, partial genome",
        "c2\tTOPOLOGY\t\tcircular\t",
        "\tsource\t1..400\torganism\tBacillus subtilis",
        "c3\tsource\t1..300\torganism\tBacillus subtilis",
        "\t\t\tff_definition\t@@[organism]@@ chromosomal gene, partial cds",
        "c4\tsource\t1..200\torganism\tBacillus subtilis",
        "\t\t\tff_definition\t@@[organism]@@ DNA, unplaced sequence @@[entry]@@",
    ]
    doc = extract(_ann_from_rows(tmp_path, rows))
    got = [(e.type, e.seq_name, e.status, e.topology) for e in doc.entries]
    assert got == [("chromosome", "", "partial", ""), ("chromosome", "", "partial", "circular"),
                   ("unplaced", "", "partial", ""), ("unplaced", "", "partial", "")]


def test_entries_follow_fasta_order(tmp_path):
    # I1: c2 has no rows (featureless contig in the COMMON-source layout); x1 is ann-only
    rows = [
        "COMMON\tsource\t1..E\torganism\tBacillus subtilis",
        "\t\t\tmol_type\tgenomic DNA",
        "c3\tCDS\t1..30\tproduct\tp3",
        "c1\tCDS\t1..30\tproduct\tp1",
        "x1\tCDS\t1..30\tproduct\tpx",
    ]
    doc = extract(_ann_from_rows(tmp_path, rows), {"c0": 10, "c1": 100, "c2": 50, "c3": 80})
    assert [e.seq_id for e in doc.entries] == ["c0", "c1", "x1", "c2", "c3"]
    # without a FASTA the ann order is kept
    assert [e.seq_id for e in extract(_ann_from_rows(tmp_path, rows)).entries] == ["c3", "c1", "x1"]


# -- per-entry DBLINK / inherited ----------------------------------------------

def test_entry_dblink_is_read(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "COMMON\tDBLINK\t\tproject\tPRJDB1",
        "A\tsource\t1..10\torganism\tX y",
        "\tDBLINK\t\tproject\tPRJDB2",
        "\t\t\tbiosample\tSAMD1",
        "B\tsource\t1..10\torganism\tX y",
    ])
    doc = extract(ann, {"A": 10, "B": 10})
    a, b = doc.entries
    assert a.dblink == {"project": "PRJDB2", "sample": "SAMD1"}
    assert b.dblink == {}
    assert doc.common["DBLINK"] == {"project": "PRJDB1"}


def test_unknown_entry_dblink_qualifier_is_noted(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "A\tsource\t1..10\torganism\tX y",
        "\tDBLINK\t\tproject\tPRJDB2",
        "\t\t\tfoo\tbar",
    ])
    doc = extract(ann, {"A": 10})
    assert Note("entry_dblink_extra", "A", "foo") in doc.notes


def test_repeated_entry_dblink_is_noted(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "A\tsource\t1..10\torganism\tX y",
        "\tDBLINK\t\tproject\tPRJDB2",
        "\tDBLINK\t\tproject\tPRJDB3",
    ])
    doc = extract(ann, {"A": 10})
    assert doc.entries[0].dblink == {"project": "PRJDB2"}
    assert Note("repeated_entry_dblink", "A", "") in doc.notes


def test_inherited_lists_common_source_names(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "COMMON\tsource\t1..E\torganism\tX y",
        "\t\t\tmol_type\tgenomic DNA",
        "\t\t\tstrain\tS0",
        "A\tsource\t1..10\tstrain\tS1",
        "\t\t\tchromosome\t1",
        "B\tTOPOLOGY\t\tcircular\t",
    ])
    doc = extract(ann, {"A": 10, "B": 10})
    a, b = doc.entries
    assert a.inherited == ["mol_type", "organism"]
    assert b.inherited == ["mol_type", "organism", "strain"]


def test_inherited_empty_without_common_source(tmp_path):
    ann = _ann_from_rows(tmp_path, ["A\tsource\t1..10\torganism\tX y"])
    assert extract(ann, {"A": 10}).entries[0].inherited == []


def test_entry_fields_round_trip():
    doc = _doc("dfast_draft")
    doc.entries[0].dblink = {"project": "PRJDB2", "sample": ["SAMD1", "SAMD2"]}
    doc.entries[0].inherited = ["organism"]
    d = doc.to_json()
    assert d["entries"][0]["dblink"] == {"project": "PRJDB2", "sample": ["SAMD1", "SAMD2"]}
    assert d["entries"][0]["inherited"] == ["organism"]
    assert MetaDoc.from_json(d) == doc


# -- where source and DBLINK are written ------------------------------------------

def test_layouts_of_fixtures():
    # name -> (source_layout, dblink_layout), decided by reading each fixture
    expected = {
        "handmade": ("common", "none"),            # COMMON source 1..E, no DBLINK
        "dfast_draft": ("entry", "common"),        # source on each entry, COMMON DBLINK
        "dfast_complete": ("entry", "common"),
        "marker_set": ("entry", "none"),
        "plasmid_only": ("entry", "common"),
        "single_chromosome": ("entry", "common"),
        "wgs_featureless": ("entry", "common"),
    }
    for name, (source_layout, dblink_layout) in expected.items():
        d = _doc(name).to_json()
        assert (d["source_layout"], d["dblink_layout"]) == (source_layout, dblink_layout), name
        assert d["layout"] in ("common_source", "per_entry")      # the old key stays
        assert (d["common_source"] != []) == (source_layout == "common"), name
        for e in d["entries"]:
            assert "source_rows" in e and "topology" in e and "dblink" in e


def test_mixed_source_layout(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "COMMON\tsource\t1..E\torganism\tX y",
        "\t\t\tmol_type\tgenomic DNA",
        "p1\tsource\t1..10\tplasmid\tpA",
        "p2\tCDS\t1..9\tproduct\thypothetical protein",
    ])
    d = extract(ann, {"p1": 10, "p2": 10}).to_json()
    assert d["source_layout"] == "mixed"
    assert d["dblink_layout"] == "none"
    assert d["common_source"] == []
    assert d["layout"] == "per_entry"
    assert [e["source_rows"] for e in d["entries"]] == [[["plasmid", "pA"]], []]


def test_partial_entry_dblink_is_mixed(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "A\tsource\t1..10\torganism\tX y",
        "\tDBLINK\t\tproject\tPRJDB2",
        "B\tsource\t1..10\torganism\tX y",
    ])
    d = extract(ann, {"A": 10, "B": 10}).to_json()
    assert d["dblink_layout"] == "mixed"
    assert d["source_layout"] == "entry"
    assert [e["dblink"] for e in d["entries"]] == [{"project": "PRJDB2"}, {}]


def test_dblink_layout_entry_common_and_both(tmp_path):
    entry_only = _ann_from_rows(tmp_path, [
        "A\tsource\t1..10\torganism\tX y",
        "\tDBLINK\t\tproject\tPRJDB2",
        "B\tsource\t1..10\torganism\tX y",
        "\tDBLINK\t\tproject\tPRJDB3",
    ])
    assert extract(entry_only, {"A": 10, "B": 10}).to_json()["dblink_layout"] == "entry"
    both = _ann_from_rows(tmp_path, [
        "COMMON\tDBLINK\t\tproject\tPRJDB1",
        "A\tsource\t1..10\torganism\tX y",
        "\tDBLINK\t\tproject\tPRJDB2",
    ])
    assert extract(both, {"A": 10}).to_json()["dblink_layout"] == "mixed"


def test_layouts_without_any_source(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "COMMON\tDBLINK\t\tproject\tPRJDB1",
        "seq1\tCDS\t1..300\tproduct\thypothetical protein",
    ])
    d = extract(ann).to_json()
    assert (d["source_layout"], d["dblink_layout"], d["common_source"]) == ("none", "common", [])


def test_common_source_raw_keeps_placeholders(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "COMMON\tsource\t1..E\torganism\tEscherichia coli",
        "\t\t\tstrain\t",
        "\t\t\tenvironmental_sample\t",
        "\t\t\tff_definition\t@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@",
        "\t\t\tsubmitter_seqid\t@@[entry]@@",
        "A\tCDS\t1..9\tproduct\thypothetical protein",
    ])
    d = extract(ann, {"A": 10}).to_json()
    assert d["source_layout"] == "common"
    assert d["common_source"] == [
        ["organism", "Escherichia coli"],
        ["strain", ""],                       # known non-flag qualifier: empty text, not true
        ["environmental_sample", True],
        ["ff_definition", "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"],
        ["submitter_seqid", "@@[entry]@@"],
    ]
    assert d["entries"][0]["source_rows"] == []


def test_entry_source_rows_in_written_order(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "A\tsource\t1..10\torganism\tX y",
        "\t\t\tnote\tfirst",
        "\t\t\tsubmitter_seqid\tA",
        "\t\t\tnote\tsecond",
        "\t\t\tff_definition\t@@[organism]@@ chromosome",
        "\t\t\tenvironmental_sample\t",
        "B\tCDS\t1..9\tproduct\thypothetical protein",
    ])
    d = extract(ann, {"A": 10, "B": 10}).to_json()
    assert d["source_layout"] == "entry"
    assert d["common_source"] == []
    assert d["entries"][0]["source_rows"] == [
        ["organism", "X y"], ["note", "first"], ["submitter_seqid", "A"], ["note", "second"],
        ["ff_definition", "@@[organism]@@ chromosome"], ["environmental_sample", True],
    ]
    assert d["entries"][1]["source_rows"] == []


def test_source_rows_round_trip():
    doc = _doc("handmade")
    doc.entries[0].source_rows = [("organism", "X y"), ("environmental_sample", True)]
    d = doc.to_json()
    assert d["entries"][0]["source_rows"] == [["organism", "X y"], ["environmental_sample", True]]
    assert MetaDoc.from_json(d) == doc
    common = _doc("handmade")
    assert common.common_source
    assert MetaDoc.from_json(common.to_json()) == common


def test_raw_rows_empty_value_is_true_only_for_flags_and_unknowns(tmp_path):
    ann = _ann_from_rows(tmp_path, [
        "COMMON\tsource\t1..E\torganism\t",
        "\t\t\tstrain\t",
        "\t\t\tenvironmental_sample\t",
        "\t\t\tmy_unknown_flag\t",
        "A\tCDS\t1..9\tproduct\thypothetical protein",
    ])
    d = extract(ann, {"A": 10}).to_json()
    assert d["common_source"] == [["organism", ""], ["strain", ""], ["environmental_sample", True],
                                  ["my_unknown_flag", True]]
    ann = _ann_from_rows(tmp_path, [
        "A\tsource\t1..10\tisolate\t",
        "\t\t\tfocus\t",
        "\t\t\tmy_unknown_flag\t",
    ])
    d = extract(ann, {"A": 10}).to_json()
    assert d["entries"][0]["source_rows"] == [["isolate", ""], ["focus", True], ["my_unknown_flag", True]]
    # the DFAST draft fixture writes empty organism / strain on each entry
    rows = dict(_doc("dfast_draft").to_json()["entries"][0]["source_rows"])
    assert rows["organism"] == "" and rows["strain"] == ""
