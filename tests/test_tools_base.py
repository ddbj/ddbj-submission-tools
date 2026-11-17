from pathlib import Path

from ddbj_gff import parse
from ddbj_gff.tools._base import ToolSpec, run_tool

FX = Path(__file__).parent / "fixtures" / "tools"
GFF = """##gff-version 3
chr1\tA\tgene\t1\t20\t.\t+\t.\tID=g1
chr1\tA\tmRNA\t1\t20\t.\t+\t.\tID=g1.t1;Parent=g1
chr1\tA\tCDS\t1\t9\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1
chr1\tA\tCDS\t15\t20\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1
chr1\tA\tintron\t10\t14\t.\t+\t.\tID=g1.t1.i1;Parent=g1.t1
chr2\tA\tgene\t1\t20\t.\t+\t.\tID=g1
chr2\tA\tmRNA\t1\t20\t.\t+\t.\tID=g1.t1;Parent=g1
chr2\tA\tCDS\t1\t9\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1
"""


def test_runs_named_passes_and_reports(tmp_path):
    src = tmp_path / "in.gff3"; src.write_text(GFF)
    out, rep = tmp_path / "out.gff3", tmp_path / "rep.txt"
    spec = ToolSpec("t", dedupe_ids=True, passes=("pass_drop_types",), config={"drop_types": ("intron",)})
    assert run_tool(spec, [str(src), "--out", str(out), "--report", str(rep)]) == 0
    text = out.read_text()
    assert "\tintron\t" not in text and "ID=chr2.g1\n" in text
    doc = parse(text)
    assert [r.id for r in doc.roots] == ["chr1.g1", "chr2.g1"]
    r = rep.read_text()
    assert "[applied]   drop-feature intron" in r and "prepass: renamed 3 id(s)" in r
    assert "##gff-version 3" in text and "##sequence-region" not in text   # no directives added here


def test_gtf_input_selects_gtf_reader(tmp_path):
    out = tmp_path / "o.gff3"
    assert run_tool(ToolSpec("t"), [str(FX / "braker3.gtf"), "--out", str(out)]) == 0
    assert "ID=g1.t1.cds" in out.read_text()


def test_snap_reader_selected_by_spec(tmp_path):
    out = tmp_path / "o.gff3"
    assert run_tool(ToolSpec("t", reader="snap"), [str(FX / "snap.gff3"), "--out", str(out)]) == 0
    assert "ID=chr1-snap.1.t1.cds" in out.read_text()


def test_missing_input_returns_2(tmp_path):
    assert run_tool(ToolSpec("t"), [str(tmp_path / "nope.gff3")]) == 2


def test_reader_diagnostics_reach_the_report(tmp_path):
    gtf = ('chr1\tA\tCDS\t1\t9\t.\t+\t0\tgene_id "g1"; transcript_id "g1.t1";\n'
           'chr1\tA\tCDS\t20\t30\t.\t+\t0\tgene_id "g2";\n'            # no transcript_id -> skipped
           'chr1\tA\tCDS\t40\t50\t.\t+\n')                               # 8 columns -> skipped
    src = tmp_path / "in.gtf"; src.write_text(gtf)
    out, rep = tmp_path / "o.gff3", tmp_path / "r.txt"
    assert run_tool(ToolSpec("t"), [str(src), "--out", str(out), "--report", str(rep)]) == 0
    r = rep.read_text()
    assert "reader: gtf-no-transcript-id 1" in r and "reader: col-count 1" in r


def test_gtf_content_is_sniffed_regardless_of_suffix(tmp_path):
    src = tmp_path / "augustus.gff"
    src.write_text((FX / "braker3.gtf").read_text())
    out = tmp_path / "o.gff3"
    assert run_tool(ToolSpec("t"), [str(src), "--out", str(out)]) == 0
    assert "ID=g1.t1.cds" in out.read_text() and "g1=" not in out.read_text()
