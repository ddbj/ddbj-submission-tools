"""nuc.mssann: .ann / FASTA reading and writing, and check.strain."""
from pathlib import Path

import pytest

from ddbj_submission_tools.nuc import mssann
from ddbj_submission_tools.nuc.check import strain as check_strain

FIX = Path(__file__).parent / "nuc_scripts_fixtures"
ANN = FIX / "small.ann"


def test_parse_reads_entries_and_common():
    common = mssann.read_common(ANN)
    assert any(f.key == "SUBMITTER" for f in common.features)
    assert [e.name for e in mssann.parse(ANN)] == ["COMMON", "seq1"]


def test_location_complement_join():
    loc = mssann.parse_location("complement(join(1..10,20..>30))")
    assert loc.strand == -1
    assert len(loc.spans) == 2


def test_writer_roundtrip(tmp_path):
    ann = mssann.read(ANN)
    out = tmp_path / "out.ann"
    with open(out, "w") as fh:
        mssann.writer.write_entries([ann.common, *ann.entries], fh)
    again = mssann.read(out)
    def shape(a):
        return [(e.name, [(f.key, f.location, [(q.key, q.value) for q in f.qualifiers]) for f in e.features])
                for e in [a.common, *a.entries]]
    assert shape(again) == shape(ann)


def test_no_host_module():
    assert not hasattr(mssann, "host")


def test_check_strain_requires_coll_dump(capsys):
    assert check_strain.main([str(ANN)]) == 2
    err = capsys.readouterr().err
    assert "coll_dump.txt" in err
    assert "https://ftp.ncbi.nih.gov/pub/taxonomy/coll_dump.txt" in err


def test_check_strain_flags_missing_space(tmp_path, capsys):
    coll = tmp_path / "coll_dump.txt"
    coll.write_text("JCM\tc\tJapan Collection of Microorganisms\nNBRC\tc\tNITE Biological Resource Center\n")
    assert check_strain.main([str(ANN), "--coll-dump", str(coll)]) == 1
    assert "'JCM 1234'" in capsys.readouterr().err
