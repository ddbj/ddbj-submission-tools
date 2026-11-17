"""`ann_meta apply --placed`: write source and DBLINK where the caller says (placed.json)."""
import json
from pathlib import Path

import pytest

from ann_fix.annfile import AnnFile
from ann_meta.cli import main
from ann_meta.extract import read_fasta_lengths

FIX = Path(__file__).parent / "ann_meta_fixtures"

# The byte-identical round trip covers the fixtures whose layout the placed JSON can describe.
# Excluded: none of the 4 is `mixed` or carries `source_location` / `repeated_source` notes
# (GGS refuses those files upstream: the JSON has no source locations / feature boundaries).
ROUNDTRIP = ["handmade", "dfast_draft", "dfast_complete", "marker_set"]
_META_KEYS = ("source", "TOPOLOGY", "DBLINK")
_COMMON_KEYS = ("SUBMITTER", "REFERENCE", "COMMENT", "ST_COMMENT", "DATE", "DBLINK", "source")


def _extract(ann: Path, fasta: Path | None, capsys) -> dict:
    argv = ["extract", "--ann", str(ann)] + (["--fasta", str(fasta)] if fasta else [])
    assert main(argv) == 0
    return json.loads(capsys.readouterr().out)


def _placed_from_doc(doc: dict) -> dict:
    """The placed JSON built from `ann_meta extract` with every layout and value as read."""
    common = doc["common"]
    out: dict = {k: common[k] for k in ("SUBMITTER", "REFERENCE", "COMMENT", "ST_COMMENT", "DATE") if k in common}
    if doc["common_extra"]:
        out["_extra"] = doc["common_extra"]
    # a file without any source is written per entry (as GGS does), every entry with no source rows
    source_layout = "entry" if doc["source_layout"] == "none" else doc["source_layout"]
    out["source_layout"] = source_layout
    out["dblink_layout"] = doc["dblink_layout"]
    if source_layout == "common":
        out["common_source"] = doc["common_source"]
    if doc["dblink_layout"] == "common":
        out["common_dblink"] = common["DBLINK"]
    entries = []
    for e in doc["entries"]:
        row = {"seq_id": e["seq_id"], "topology": e["topology"]}
        if source_layout == "entry":
            row["source"] = e["source_rows"]
        if doc["dblink_layout"] == "entry":
            row["dblink"] = e["dblink"]
        entries.append(row)
    out["entries"] = entries
    return out


def _apply(tmp_path: Path, ann: Path, placed: dict, *, fasta: Path | None = None, out: Path | None = None,
           extra_args: list[str] | None = None) -> tuple[int, Path]:
    pj = tmp_path / "placed.json"
    pj.write_text(json.dumps(placed, ensure_ascii=False), encoding="utf-8")
    out = out or tmp_path / "out.ann"
    argv = ["apply", "--placed", "--ann", str(ann), "--common", str(pj), "--out", str(out)]
    if fasta:
        argv += ["--fasta", str(fasta)]
    return main(argv + (extra_args or [])), out


def _feature_rows(path: Path) -> list[tuple[str, bytes]]:
    """(entry, line without column 1) of every row outside the rows apply --placed may rewrite."""
    out = []
    entry = key = ""
    for raw in path.read_bytes().splitlines(keepends=True):
        cols = raw.split(b"\t", 1)
        if cols[0]:
            entry = cols[0].decode()
        rest = cols[1] if len(cols) > 1 else b""
        k = rest.split(b"\t", 1)[0].rstrip(b"\r\n").decode()
        if k:
            key = k
        if key in _META_KEYS or (entry == "COMMON" and key in _COMMON_KEYS):
            continue
        out.append((entry, rest))
    return out


def _features(path: Path):
    ann = AnnFile.read(path)
    return ann, ann.features()


# -- round trip ------------------------------------------------------------------

def _variants(data: bytes):
    lf = data.replace(b"\r\n", b"\n")
    crlf = lf.replace(b"\n", b"\r\n")
    return {"as_is": data, "lf": lf, "crlf": crlf, "lf_nofinal": lf.rstrip(b"\n"),
            "crlf_nofinal": crlf.rstrip(b"\r\n")}


@pytest.mark.parametrize("name", ROUNDTRIP)
def test_placed_roundtrip_byte_identical(name, tmp_path, capsys):
    fasta = FIX / f"{name}.fasta"
    for variant, body in _variants((FIX / f"{name}.ann").read_bytes()).items():
        src = tmp_path / f"{variant}.ann"
        src.write_bytes(body)
        doc = _extract(src, fasta, capsys)
        assert doc["source_layout"] in ("common", "entry") and doc["dblink_layout"] != "mixed"
        assert not {n["code"] for n in doc["notes"]} & {"source_location", "repeated_source"}
        rc, out = _apply(tmp_path, src, _placed_from_doc(doc), fasta=fasta, out=tmp_path / f"{variant}.out.ann")
        assert rc == 0, capsys.readouterr()
        assert out.read_bytes() == body, (name, variant)


def test_placed_roundtrip_without_fasta(tmp_path, capsys):
    src = FIX / "dfast_complete.ann"
    doc = _extract(src, None, capsys)
    rc, out = _apply(tmp_path, src, _placed_from_doc(doc))
    assert rc == 0
    assert out.read_bytes() == src.read_bytes()


# -- source layout -----------------------------------------------------------------

def test_placed_common_to_entry_source(tmp_path, capsys):
    src, fasta = FIX / "handmade.ann", FIX / "handmade.fasta"
    lengths = read_fasta_lengths(fasta)
    doc = _extract(src, fasta, capsys)
    assert doc["source_layout"] == "common"
    placed = _placed_from_doc(doc)
    placed["source_layout"] = "entry"
    del placed["common_source"]
    for e in placed["entries"]:
        e["source"] = [["organism", "Bacillus subtilis"], ["mol_type", "genomic DNA"], ["strain", "HS-7"],
                       ["environmental_sample", True], ["note", "a"], ["note", "b"],
                       ["ff_definition", f"@@[organism]@@ @@[strain]@@ DNA, contig: {e['seq_id']}"]]
    rc, out = _apply(tmp_path, src, placed, fasta=fasta)
    assert rc == 0, capsys.readouterr()
    ann, feats = _features(out)
    assert not [f for f in feats if f.entry == "COMMON" and f.key == "source"]
    sources = {f.entry: f for f in feats if f.key == "source"}
    assert set(sources) == set(lengths)
    for sid, f in sources.items():
        assert f.location == f"1..{lengths[sid]}"
        assert f.qualifiers(ann.rows) == [
            ("organism", "Bacillus subtilis"), ("mol_type", "genomic DNA"), ("strain", "HS-7"),
            ("environmental_sample", ""), ("note", "a"), ("note", "b"),
            ("ff_definition", f"@@[organism]@@ @@[strain]@@ DNA, contig: {sid}")]
        # the source is the first feature of its entry and carries the entry marker
        assert ann.rows[f.first].cols[0] == sid
    assert _feature_rows(out) == _feature_rows(src)
    data = out.read_bytes()
    assert data.count(b"\n") == data.count(b"\r\n")          # CRLF file: no bare LF
    assert b"\t\t\tenvironmental_sample\t\r\n" in data       # valueless qualifier


def test_placed_entry_to_common_source(tmp_path, capsys):
    src, fasta = FIX / "dfast_draft.ann", FIX / "dfast_draft.fasta"
    doc = _extract(src, fasta, capsys)
    assert doc["source_layout"] == "entry"
    placed = _placed_from_doc(doc)
    placed["source_layout"] = "common"
    placed["common_source"] = [["organism", "Escherichia coli"], ["mol_type", "genomic DNA"],
                               ["submitter_seqid", "@@[entry]@@"],
                               ["ff_definition", "@@[organism]@@ DNA, @@[submitter_seqid]@@"]]
    for e in placed["entries"]:
        e.pop("source", None)
    rc, out = _apply(tmp_path, src, placed, fasta=fasta)
    assert rc == 0, capsys.readouterr()
    ann, feats = _features(out)
    common_src = [f for f in feats if f.entry == "COMMON" and f.key == "source"]
    assert len(common_src) == 1 and common_src[0].location == "1..E"
    assert common_src[0].qualifiers(ann.rows) == [tuple(p) for p in placed["common_source"]]
    assert not [f for f in feats if f.entry != "COMMON" and f.key == "source"]
    assert _feature_rows(out) == _feature_rows(src)
    # each entry keeps its marker on its first remaining row
    assert [r.cols[0] for r in ann.rows if r.cols[0]] == ["COMMON", "sequence1", "sequence2", "sequence3"]
    assert b"\r" not in out.read_bytes()


def test_placed_entry_topology_and_new_entry_from_fasta(tmp_path, capsys):
    """A FASTA entry without rows gets its source (and TOPOLOGY) at its FASTA position."""
    src = tmp_path / "in.ann"
    src.write_bytes(b"COMMON\tSUBMITTER\t\tcontact\tA B\n"
                    b"s1\tsource\t1..10\torganism\tX\n"
                    b"s3\tsource\t1..30\torganism\tX\n"
                    b"\tCDS\t1..30\tproduct\tp")
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">s1\n" + "A" * 10 + "\n>s2\n" + "A" * 20 + "\n>s3\n" + "A" * 30 + "\n")
    placed = {"SUBMITTER": {"contact": "A B"}, "source_layout": "entry", "dblink_layout": "none",
              "entries": [{"seq_id": s, "topology": t, "source": [["organism", "X"]]}
                          for s, t in (("s1", ""), ("s2", "circular"), ("s3", "circular"))]}
    rc, out = _apply(tmp_path, src, placed, fasta=fasta)
    assert rc == 0, capsys.readouterr()
    assert out.read_bytes() == (b"COMMON\tSUBMITTER\t\tcontact\tA B\n"
                                b"s1\tsource\t1..10\torganism\tX\n"
                                b"s2\tTOPOLOGY\t\tcircular\t\n"
                                b"\tsource\t1..20\torganism\tX\n"
                                b"s3\tTOPOLOGY\t\tcircular\t\n"
                                b"\tsource\t1..30\torganism\tX\n"
                                b"\tCDS\t1..30\tproduct\tp")


def test_placed_common_layout_topology_and_removed_topology(tmp_path, capsys):
    src, fasta = FIX / "dfast_complete.ann", FIX / "dfast_complete.fasta"
    doc = _extract(src, fasta, capsys)
    placed = _placed_from_doc(doc)
    placed["source_layout"] = "common"
    placed["common_source"] = [["organism", "Escherichia coli"], ["mol_type", "genomic DNA"]]
    placed["entries"][0]["topology"] = "circular"
    placed["entries"][1]["topology"] = ""
    for e in placed["entries"]:
        e.pop("source", None)
    rc, out = _apply(tmp_path, src, placed, fasta=fasta)
    assert rc == 0, capsys.readouterr()
    ann, feats = _features(out)
    topo = [f.entry for f in feats if f.key == "TOPOLOGY"]
    assert topo == ["sequence1"]
    first = next(f for f in feats if f.entry == "sequence1")
    assert first.key == "TOPOLOGY" and ann.rows[first.first].cols[0] == "sequence1"
    assert _feature_rows(out) == _feature_rows(src)


# -- DBLINK ---------------------------------------------------------------------------

def test_placed_dblink_none_removes_both(tmp_path, capsys):
    src = tmp_path / "in.ann"
    src.write_bytes(b"COMMON\tDBLINK\t\tproject\tPRJDB1\n"
                    b"\tSUBMITTER\t\tcontact\tA B\n"
                    b"s1\tDBLINK\t\tbiosample\tSAMD1\n"
                    b"\tsource\t1..10\torganism\tX\n"
                    b"\tCDS\t1..9\tproduct\tp\n")
    placed = {"SUBMITTER": {"contact": "A B"}, "source_layout": "entry", "dblink_layout": "none",
              "entries": [{"seq_id": "s1", "topology": "", "source": [["organism", "X"]]}]}
    rc, out = _apply(tmp_path, src, placed)
    assert rc == 0, capsys.readouterr()
    assert out.read_bytes() == (b"COMMON\tSUBMITTER\t\tcontact\tA B\n"
                                b"s1\tsource\t1..10\torganism\tX\n"
                                b"\tCDS\t1..9\tproduct\tp\n")


def test_placed_dblink_entry_writes_before_source(tmp_path, capsys):
    src, fasta = FIX / "dfast_complete.ann", FIX / "dfast_complete.fasta"
    doc = _extract(src, fasta, capsys)
    assert doc["dblink_layout"] == "common"
    placed = _placed_from_doc(doc)
    placed["dblink_layout"] = "entry"
    del placed["common_dblink"]
    placed["entries"][0]["dblink"] = {"project": "PRJDB99901", "sample": ["SAMD99000001"],
                                      "sequence_read_archive": ["DRR990001", "DRR990002"]}
    placed["entries"][1]["dblink"] = {"project": "PRJDB99901", "sample": "SAMD99000002"}
    rc, out = _apply(tmp_path, src, placed, fasta=fasta)
    assert rc == 0, capsys.readouterr()
    ann, feats = _features(out)
    assert not [f for f in feats if f.entry == "COMMON" and f.key == "DBLINK"]
    for sid, quals in (("sequence1", [("project", "PRJDB99901"), ("biosample", "SAMD99000001"),
                                      ("sequence read archive", "DRR990001"),
                                      ("sequence read archive", "DRR990002")]),
                       ("sequence2", [("project", "PRJDB99901"), ("biosample", "SAMD99000002")])):
        ef = [f for f in feats if f.entry == sid]
        i = next(i for i, f in enumerate(ef) if f.key == "DBLINK")
        assert ef[i + 1].key == "source"
        assert ef[i].qualifiers(ann.rows) == quals
        assert ann.rows[ef[0].first].cols[0] == sid
    assert _feature_rows(out) == _feature_rows(src)
    # the other COMMON rows are unchanged
    assert out.read_bytes().startswith(b"COMMON\tSUBMITTER\t\tab_name\tMishima,T.\n")


def test_placed_dblink_common_from_entry(tmp_path, capsys):
    src = tmp_path / "in.ann"
    src.write_bytes(b"COMMON\tSUBMITTER\t\tcontact\tA B\n"
                    b"s1\tDBLINK\t\tbiosample\tSAMD1\n"
                    b"\tsource\t1..10\torganism\tX\n")
    placed = {"SUBMITTER": {"contact": "A B"}, "source_layout": "entry", "dblink_layout": "common",
              "common_dblink": {"project": "PRJDB1", "sample": "SAMD1"},
              "entries": [{"seq_id": "s1", "topology": "", "source": [["organism", "X"]]}]}
    rc, out = _apply(tmp_path, src, placed)
    assert rc == 0, capsys.readouterr()
    assert out.read_bytes() == (b"COMMON\tDBLINK\t\tproject\tPRJDB1\n"
                                b"\t\t\tbiosample\tSAMD1\n"
                                b"\tSUBMITTER\t\tcontact\tA B\n"
                                b"s1\tsource\t1..10\torganism\tX\n")


# -- COMMON ------------------------------------------------------------------------------

def test_placed_keeps_date(tmp_path, capsys):
    src, fasta = FIX / "handmade.ann", FIX / "handmade.fasta"
    doc = _extract(src, fasta, capsys)
    placed = _placed_from_doc(doc)
    placed["DATE"] = {"hold_date": "20991231"}
    placed["_extra"] = {**placed.get("_extra", {}), "DATE": [["note", "x"]]}
    placed["REFERENCE"][0]["title"] = "A new title"
    placed["COMMENT"] = [{"line": ["first", "second"]}]
    rc, out = _apply(tmp_path, src, placed, fasta=fasta)
    assert rc == 0, capsys.readouterr()
    data = out.read_bytes()
    assert data.count(b"\tDATE\t") == 1
    assert b"\tDATE\t\thold_date\t20271231\r\n" in data
    assert b"20991231" not in data and b"\tnote\tx" not in data
    ann, feats = _features(out)
    keys = [f.key for f in feats if f.entry == "COMMON"]
    assert keys == ["SUBMITTER", "REFERENCE", "COMMENT", "DATE", "KEYWORD", "source"]
    ref = next(f for f in feats if f.key == "REFERENCE")
    assert ref.get(ann.rows, "title") == "A new title"
    sub = next(f for f in feats if f.key == "SUBMITTER")
    assert sub.get(ann.rows, "fax") == "81-55-981-0000"


def test_placed_removes_comment_and_st_comment_when_absent(tmp_path, capsys):
    src, fasta = FIX / "dfast_complete.ann", FIX / "dfast_complete.fasta"
    doc = _extract(src, fasta, capsys)
    placed = _placed_from_doc(doc)
    del placed["COMMENT"]
    placed["ST_COMMENT"] = []
    rc, out = _apply(tmp_path, src, placed, fasta=fasta)
    assert rc == 0, capsys.readouterr()
    keys = [f.key for f in AnnFile.read(out).features() if f.entry == "COMMON"]
    assert keys == ["DBLINK", "SUBMITTER", "REFERENCE"]
    assert _feature_rows(out) == _feature_rows(src)


def _roundtrip_body(tmp_path, capsys, body: bytes, mutate=None) -> bytes:
    src = tmp_path / "in.ann"
    src.write_bytes(body)
    placed = _placed_from_doc(_extract(src, None, capsys))
    if mutate:
        mutate(placed)
    rc, out = _apply(tmp_path, src, placed)
    assert rc == 0, capsys.readouterr()
    return out.read_bytes()


@pytest.mark.parametrize("body", [
    # SUBMITTER before DBLINK
    b"COMMON\tSUBMITTER\t\tcontact\tA B\n\t\t\temail\ta@b.c\n\tDBLINK\t\tproject\tPRJDB1\n"
    b"s1\tsource\t1..10\torganism\tX\n\tCDS\t1..9\tproduct\tp\n",
    # ST_COMMENT before COMMENT, DATE first
    b"COMMON\tDATE\t\thold_date\t20271231\n\tSUBMITTER\t\tcontact\tA B\n"
    b"\tST_COMMENT\t\ttagset_id\tGenome-Assembly-Data\n\t\t\tAssembly Method\tX\n\tCOMMENT\t\tline\thello\n"
    b"s1\tsource\t1..10\torganism\tX\n",
    # REFERENCE before SUBMITTER, two COMMENTs around them
    b"COMMON\tCOMMENT\t\tline\ta\n\tREFERENCE\t\ttitle\tT\n\t\t\tab_name\tA,B.\n\t\t\tstatus\tUnpublished\n"
    b"\tSUBMITTER\t\tcontact\tA B\n\tCOMMENT\t\tline\tb\ns1\tsource\t1..10\torganism\tX\n",
])
def test_placed_noncanonical_common_order_noop(tmp_path, capsys, body):
    assert _roundtrip_body(tmp_path, capsys, body) == body


def _rows_of(data: bytes, entry: str, keys: tuple[str, ...]) -> list[bytes]:
    """Lines (without column 1) of the *keys* features of *entry*."""
    out, cur, key = [], "", ""
    for raw in data.splitlines(keepends=True):
        cols = raw.split(b"\t")
        if cols[0]:
            cur = cols[0].decode()
        if len(cols) > 1 and cols[1]:
            key = cols[1].decode()
        if cur == entry and key in keys:
            out.append(raw.split(b"\t", 1)[1])
    return out


def test_placed_dblink_switch_keeps_other_common_rows(tmp_path, capsys):
    """Switching the DBLINK layout or dropping COMMENT leaves SUBMITTER / REFERENCE / ST_COMMENT as they were
    (fax not moved, DFAST's empty placeholders kept)."""
    src = FIX / "dfast_draft.ann"
    body = src.read_bytes()

    def mutate(p):
        p["dblink_layout"] = "none"
        p.pop("common_dblink", None)
        p.pop("COMMENT", None)
    out = _roundtrip_body(tmp_path, capsys, body, mutate)
    keys = ("SUBMITTER", "REFERENCE", "ST_COMMENT")
    assert _rows_of(out, "COMMON", keys) == _rows_of(body, "COMMON", keys)
    assert b"\tDBLINK\t" not in out and b"\tCOMMENT\t" not in out

    fax = (b"COMMON\tDBLINK\t\tproject\tPRJDB1\n\tSUBMITTER\t\tcontact\tA B\n\t\t\tfax\t1\n\t\t\temail\ta@b.c\n"
           b"s1\tsource\t1..10\torganism\tX\n")
    out = _roundtrip_body(tmp_path, capsys, fax, mutate)
    assert out == (b"COMMON\tSUBMITTER\t\tcontact\tA B\n\t\t\tfax\t1\n\t\t\temail\ta@b.c\n"
                   b"s1\tsource\t1..10\torganism\tX\n")


def test_placed_new_common_feature_goes_to_canonical_place(tmp_path, capsys):
    body = (b"COMMON\tSUBMITTER\t\tcontact\tA B\n\tREFERENCE\t\ttitle\tT\n"
            b"\tST_COMMENT\t\ttagset_id\tGenome-Assembly-Data\n\t\t\tAssembly Method\tX\n"
            b"s1\tsource\t1..10\torganism\tX\n")

    def mutate(p):
        p["COMMENT"] = [{"line": "new"}]
        p["dblink_layout"] = "common"
        p["common_dblink"] = {"project": "PRJDB1"}
    out = _roundtrip_body(tmp_path, capsys, body, mutate)
    assert out == (b"COMMON\tDBLINK\t\tproject\tPRJDB1\n\tSUBMITTER\t\tcontact\tA B\n\tREFERENCE\t\ttitle\tT\n"
                   b"\tCOMMENT\t\tline\tnew\n"
                   b"\tST_COMMENT\t\ttagset_id\tGenome-Assembly-Data\n\t\t\tAssembly Method\tX\n"
                   b"s1\tsource\t1..10\torganism\tX\n")


def test_placed_changed_common_feature_rewritten_in_place(tmp_path, capsys):
    body = (b"COMMON\tREFERENCE\t\ttitle\tT\n\tSUBMITTER\t\tcontact\tA B\n\t\t\tfax\t1\n\t\t\temail\ta@b.c\n"
            b"s1\tsource\t1..10\torganism\tX\n")
    out = _roundtrip_body(tmp_path, capsys, body, lambda p: p["REFERENCE"][0].update(title="U"))
    assert out == body.replace(b"title\tT", b"title\tU")


def test_placed_keeps_linear_topology(tmp_path, capsys):
    body = b"COMMON\tSUBMITTER\t\tcontact\tA B\ns1\tTOPOLOGY\t\tlinear\t\n\tsource\t1..10\torganism\tX\n"
    assert _roundtrip_body(tmp_path, capsys, body) == body
    # circular asked for: the linear row is replaced
    out = _roundtrip_body(tmp_path, capsys, body, lambda p: p["entries"][0].update(topology="circular"))
    assert out == body.replace(b"linear", b"circular")


_UNSORTED_ST = (b"COMMON\tSUBMITTER\t\tcontact\tA B\n\t\t\tfax\t1\n\t\t\temail\ta@b.c\n"
                b"\tST_COMMENT\t\ttagset_id\tGenome-Assembly-Data\n\t\t\tGenome Coverage\t10x\n"
                b"\t\t\tAssembly Method\tX\n\t\t\tSequencing Technology\tIllumina\n"
                b"s1\tsource\t1..10\torganism\tX\n")


def test_placed_unsorted_st_comment_noop_is_byte_identical(tmp_path, capsys):
    assert _roundtrip_body(tmp_path, capsys, _UNSORTED_ST) == _UNSORTED_ST


def test_placed_st_comment_edit_keeps_the_input_tag_order(tmp_path, capsys):
    """A real ST_COMMENT edit writes the tags in the input's order (not create_common's sorted order)."""
    def edit(p):
        p["ST_COMMENT"]["Genome Coverage"] = "20x"
    out = _roundtrip_body(tmp_path, capsys, _UNSORTED_ST, edit)
    assert out == _UNSORTED_ST.replace(b"Genome Coverage\t10x", b"Genome Coverage\t20x")

    def reorder(p):                      # a tag order change is a real change: written as given
        st = p["ST_COMMENT"]
        p["ST_COMMENT"] = {"tagset_id": st["tagset_id"], "Assembly Method": st["Assembly Method"],
                           "Sequencing Technology": st["Sequencing Technology"],
                           "Genome Coverage": st["Genome Coverage"]}
    out = _roundtrip_body(tmp_path, capsys, _UNSORTED_ST, reorder)
    assert out == _UNSORTED_ST.replace(
        b"\t\t\tGenome Coverage\t10x\n\t\t\tAssembly Method\tX\n\t\t\tSequencing Technology\tIllumina\n",
        b"\t\t\tAssembly Method\tX\n\t\t\tSequencing Technology\tIllumina\n\t\t\tGenome Coverage\t10x\n")


def test_placed_st_comment_list_values_and_two_sections(tmp_path, capsys):
    body = (b"COMMON\tSUBMITTER\t\tcontact\tA B\n"
            b"\tST_COMMENT\t\ttagset_id\tA\n\t\t\tz\t1\n"
            b"\tST_COMMENT\t\ttagset_id\tB\n\t\t\ty\t2\n"
            b"s1\tsource\t1..10\torganism\tX\n")

    def edit(p):
        p["ST_COMMENT"][1]["y"] = ["2", "3"]
    out = _roundtrip_body(tmp_path, capsys, body, edit)
    assert out == body.replace(b"\ty\t2\n", b"\ty\t2; 3\n")


def test_placed_keeps_placeholders_and_ff_definition(tmp_path, capsys):
    src, fasta = FIX / "plasmid_only.ann", FIX / "plasmid_only.fasta"
    doc = _extract(src, fasta, capsys)
    placed = _placed_from_doc(doc)
    ff = "@@[organism]@@ @@[strain]@@ plasmid @@[plasmid]@@ DNA, complete sequence"
    placed["entries"][0]["source"] = [["organism", "Escherichia coli"], ["strain", "NIG-0002b"],
                                      ["plasmid", "pNIG1"], ["ff_definition", ff]]
    rc, out = _apply(tmp_path, src, placed, fasta=fasta)
    assert rc == 0, capsys.readouterr()
    ann, feats = _features(out)
    s = next(f for f in feats if f.entry == "pA" and f.key == "source")
    # written exactly as given: no mol_type added, no ff_definition regenerated, placeholders kept
    assert s.qualifiers(ann.rows) == [("organism", "Escherichia coli"), ("strain", "NIG-0002b"),
                                      ("plasmid", "pNIG1"), ("ff_definition", ff)]
    assert _feature_rows(out) == _feature_rows(src)


# -- entries without a source -----------------------------------------------------------

_PARTIAL_SOURCE = (b"COMMON\tSUBMITTER\t\tcontact\tA B\n"
                   b"c1\tsource\t1..10\torganism\tX\n"
                   b"\t\t\tmol_type\tgenomic DNA\n"
                   b"\tCDS\t1..9\tproduct\tp\n"
                   b"c2\tTOPOLOGY\t\tcircular\t\n"
                   b"\tCDS\t1..20\tproduct\tq\n")
_NO_SOURCE = (b"COMMON\tSUBMITTER\t\tcontact\tA B\n"
              b"c1\tCDS\t1..9\tproduct\tp\n"
              b"c2\tTOPOLOGY\t\tcircular\t\n"
              b"\tCDS\t1..20\tproduct\tq\n")


def _fasta(tmp_path: Path, *entries: tuple[str, int]) -> Path:
    fasta = tmp_path / "in.fasta"
    fasta.write_text("".join(f">{sid}\n{'A' * n}\n" for sid, n in entries))
    return fasta


def _roundtrip_fasta(tmp_path, capsys, body: bytes, fasta: Path, mutate=None) -> tuple[dict, bytes]:
    src = tmp_path / "in.ann"
    src.write_bytes(body)
    placed = _placed_from_doc(_extract(src, fasta, capsys))
    if mutate:
        mutate(placed)
    rc, out = _apply(tmp_path, src, placed, fasta=fasta)
    assert rc == 0, capsys.readouterr()
    return placed, out.read_bytes()


def test_placed_entry_without_source_noop_is_byte_identical(tmp_path, capsys):
    """c1 has a source, c2 none: c2 comes back with no source rows and stays without a source."""
    fasta = _fasta(tmp_path, ("c1", 10), ("c2", 20))
    for variant, body in _variants(_PARTIAL_SOURCE).items():
        placed, out = _roundtrip_fasta(tmp_path, capsys, body, fasta)
        assert placed["source_layout"] == "entry"
        assert [e["source"] for e in placed["entries"]] == [[["organism", "X"], ["mol_type", "genomic DNA"]], []]
        assert out == body, variant


def test_placed_file_without_source_noop_is_byte_identical(tmp_path, capsys):
    """A file without any source: layout entry with no source rows for every entry changes nothing."""
    fasta = _fasta(tmp_path, ("c1", 10), ("c2", 20))
    for variant, body in _variants(_NO_SOURCE).items():
        placed, out = _roundtrip_fasta(tmp_path, capsys, body, fasta)
        assert placed["source_layout"] == "entry"
        assert placed["entries"] == [{"seq_id": "c1", "topology": "", "source": []},
                                     {"seq_id": "c2", "topology": "circular", "source": []}]
        assert out == body, variant


def test_placed_adds_a_source_to_an_entry_without_one(tmp_path, capsys):
    fasta = _fasta(tmp_path, ("c1", 10), ("c2", 20))
    _, out = _roundtrip_fasta(tmp_path, capsys, _PARTIAL_SOURCE, fasta,
                              lambda p: p["entries"][1].update(source=[["organism", "Y"]]))
    assert out == _PARTIAL_SOURCE.replace(b"c2\tTOPOLOGY\t\tcircular\t\n",
                                          b"c2\tTOPOLOGY\t\tcircular\t\n\tsource\t1..20\torganism\tY\n")
    _, out = _roundtrip_fasta(tmp_path, capsys, _NO_SOURCE, fasta,
                              lambda p: p["entries"][0].update(source=[["organism", "Y"]]))
    assert out == _NO_SOURCE.replace(b"c1\tCDS\t1..9\tproduct\tp\n",
                                     b"c1\tsource\t1..10\torganism\tY\n\tCDS\t1..9\tproduct\tp\n")


def test_placed_empty_source_removes_an_existing_one(tmp_path, capsys):
    """An entry that had a source and is given no source rows loses it (an explicit edit)."""
    fasta = _fasta(tmp_path, ("c1", 10), ("c2", 20))
    _, out = _roundtrip_fasta(tmp_path, capsys, _PARTIAL_SOURCE, fasta,
                              lambda p: p["entries"][0].update(source=[]))
    assert out == _PARTIAL_SOURCE.replace(b"c1\tsource\t1..10\torganism\tX\n\t\t\tmol_type\tgenomic DNA\n"
                                          b"\tCDS", b"c1\tCDS")


def test_placed_common_source_still_required(tmp_path, capsys):
    placed = _placed(tmp_path, capsys, "handmade")
    assert placed["source_layout"] == "common"
    placed["common_source"] = []
    issues = _rejected(tmp_path, capsys, placed, "handmade")
    assert {"field": "common_source", "entry": None, "message": "missing_source"} in issues


_DUP_ST = (b"COMMON\tSUBMITTER\t\tcontact\tA B\n"
           b"\tST_COMMENT\t\ttagset_id\tGenome-Assembly-Data\n"
           b"\t\t\tSequencing Technology\tIllumina\n"
           b"\t\t\tSequencing Technology\tNanopore\n"
           b"\t\t\tAssembly Method\tFlye\n"
           b"s1\tsource\t1..10\torganism\tX\n")


def test_placed_repeated_st_comment_tag_noop_is_byte_identical(tmp_path, capsys):
    """A repeated ST_COMMENT tag comes back in _extra; it is written next to the tag it repeats."""
    src = tmp_path / "in.ann"
    src.write_bytes(_DUP_ST)
    doc = _extract(src, None, capsys)
    assert doc["common_extra"] == {"ST_COMMENT": [["Sequencing Technology", "Nanopore"]]}
    assert _roundtrip_body(tmp_path, capsys, _DUP_ST) == _DUP_ST

    def edit(p):
        p["ST_COMMENT"]["Assembly Method"] = "Canu"
    assert _roundtrip_body(tmp_path, capsys, _DUP_ST, edit) == _DUP_ST.replace(b"Flye", b"Canu")


# -- rejections ----------------------------------------------------------------------------

def _rejected(tmp_path, capsys, placed, name="dfast_complete") -> list[dict]:
    capsys.readouterr()
    out = tmp_path / "rejected.ann"
    rc, _ = _apply(tmp_path, FIX / f"{name}.ann", placed, fasta=FIX / f"{name}.fasta", out=out)
    assert rc == 2
    assert not out.exists()
    return json.loads(capsys.readouterr().out)["issues"]


def _placed(tmp_path, capsys, name="dfast_complete") -> dict:
    return _placed_from_doc(_extract(FIX / f"{name}.ann", FIX / f"{name}.fasta", capsys))


def test_placed_rejects_unknown_entry(tmp_path, capsys):
    placed = _placed(tmp_path, capsys)
    placed["entries"].append({"seq_id": "nope", "topology": "", "source": [["organism", "X"]]})
    issues = _rejected(tmp_path, capsys, placed)
    assert {"field": "entries", "entry": "nope", "message": "unknown_entry"} in issues


def test_placed_rejects_missing_entry(tmp_path, capsys):
    placed = _placed(tmp_path, capsys)
    placed["entries"].pop()
    issues = _rejected(tmp_path, capsys, placed)
    assert {"field": "entries", "entry": "sequence2", "message": "missing_entry"} in issues


def test_placed_rejects_duplicate_entry(tmp_path, capsys):
    placed = _placed(tmp_path, capsys)
    placed["entries"].append(dict(placed["entries"][0]))
    issues = _rejected(tmp_path, capsys, placed)
    assert any(i["message"] == "duplicate_seq_id" and i["entry"] == "sequence1" for i in issues)


@pytest.mark.parametrize("bad", ["E. coli\x00", "a\nb", "a\tb", "a b", "a\x7fb"])
def test_placed_rejects_bad_chars(tmp_path, capsys, bad):
    placed = _placed(tmp_path, capsys)
    placed["entries"][0]["source"][0] = ["organism", bad]
    issues = _rejected(tmp_path, capsys, placed)
    assert any(i["message"] == "forbidden_char" for i in issues)


def test_placed_rejects_bad_chars_in_common(tmp_path, capsys):
    placed = _placed(tmp_path, capsys)
    placed["SUBMITTER"]["institute"] = "x\ry"
    issues = _rejected(tmp_path, capsys, placed)
    assert {"field": "SUBMITTER.institute", "entry": None, "message": "forbidden_char"} in issues


@pytest.mark.parametrize("mutate, message", [
    (lambda p: p.update(source_layout="mixed"), "bad_source_layout"),
    (lambda p: p.update(dblink_layout="both"), "bad_dblink_layout"),
    (lambda p: p.pop("source_layout"), "bad_source_layout"),
    (lambda p: p["entries"][0].update(topology="linear"), "bad_topology"),
    (lambda p: p["entries"][0]["source"].append(["bad name", "x"]), "bad_qualifier_name"),
    (lambda p: p["entries"][0]["source"].append(["note", False]), "bad_shape"),
    (lambda p: p["entries"][0]["source"].append(["note"]), "bad_shape"),
    (lambda p: p["entries"][0].update(source="organism"), "bad_shape"),
    (lambda p: p["entries"][0].update(extra=1), "not_allowed_with_placed"),
    (lambda p: p.update(KEYWORD={"keyword": "WGS"}), "not_allowed_with_placed"),
    (lambda p: p.update(SUBMITTER={"ab_name": "not a list"}), "bad_shape"),
    (lambda p: p.update(common_dblink={"project": ["a", "b"]}), "bad_shape"),
    (lambda p: p.update(entries={}), "bad_shape"),
    (lambda p: p["entries"][0].update(seq_id=""), "empty_seq_id"),
])
def test_placed_rejects_bad_shape(tmp_path, capsys, mutate, message):
    placed = _placed(tmp_path, capsys)
    mutate(placed)
    issues = _rejected(tmp_path, capsys, placed)
    assert any(i["message"].startswith(message) for i in issues), issues


def test_placed_rejects_new_source_without_length(tmp_path, capsys):
    """Without --fasta the length comes from the original 1..N; an entry without one cannot get a source."""
    src = tmp_path / "in.ann"
    src.write_bytes(b"COMMON\tsource\t1..E\torganism\tX\n"
                    b"s1\tCDS\t1..9\tproduct\tp\n")
    placed = {"source_layout": "entry", "dblink_layout": "none",
              "entries": [{"seq_id": "s1", "topology": "", "source": [["organism", "X"]]}]}
    capsys.readouterr()
    rc, out = _apply(tmp_path, src, placed)
    assert rc == 2 and not out.exists()
    issues = json.loads(capsys.readouterr().out)["issues"]
    assert {"field": "entries", "entry": "s1", "message": "length_unknown"} in issues


def test_placed_flags_are_exclusive(tmp_path, capsys):
    placed = _placed(tmp_path, capsys)
    src = FIX / "dfast_complete.ann"
    for extra in (["--entries", str(tmp_path / "e.tsv")], ["--category", "GNM"], ["--common-only"],
                  ["--allow-unknown", "x"]):
        out = tmp_path / "x.ann"
        rc, _ = _apply(tmp_path, src, placed, out=out, extra_args=extra)
        assert rc == 1, extra
        assert not out.exists()
        assert capsys.readouterr().out == ""


def test_placed_out_may_equal_ann(tmp_path, capsys):
    src = tmp_path / "in.ann"
    src.write_bytes((FIX / "dfast_complete.ann").read_bytes())
    placed = _placed(tmp_path, capsys)
    placed["REFERENCE"][0]["title"] = "Changed"
    rc, out = _apply(tmp_path, src, placed, fasta=FIX / "dfast_complete.fasta", out=src)
    assert rc == 0
    ann, feats = _features(src)
    assert next(f for f in feats if f.key == "REFERENCE").get(ann.rows, "title") == "Changed"


def test_placed_drops_fax_when_extra_no_longer_has_it(tmp_path, capsys):
    """GGS drops SUBMITTER phone / fax (not in the spec) from _extra; the rewritten SUBMITTER then has no fax row."""
    body = (b"COMMON\tSUBMITTER\t\tcontact\tA B\n\t\t\tfax\t1\n\t\t\temail\ta@b.c\n"
            b"s1\tsource\t1..10\torganism\tX\n")

    def mutate(p):
        p["_extra"].pop("SUBMITTER", None)
    out = _roundtrip_body(tmp_path, capsys, body, mutate)
    assert out == (b"COMMON\tSUBMITTER\t\tcontact\tA B\n\t\t\temail\ta@b.c\n"
                   b"s1\tsource\t1..10\torganism\tX\n")
