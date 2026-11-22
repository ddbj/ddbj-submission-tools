"""nuc.check.organelle_ref end to end, with synthetic GenBank references (no network)."""
import json
import random

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from ddbj_submission_tools.nuc.check import organelle_ref as cli
from ddbj_submission_tools.nuc.fix import edit
from ddbj_submission_tools.nuc.organelle_ref import align, refs

needs_blast = pytest.mark.skipif(not align.blast_available(), reason="blastn not on PATH")

_rng = random.Random(11)
CERANA = "".join(_rng.choice("ACGT") for _ in range(3000))
MELLIFERA = "".join(("A" if b != "A" else "C") if i % 20 == 0 else b for i, b in enumerate(CERANA))   # 95%


def write_ref(path, acc, species, seq):
    rec = SeqRecord(Seq(seq), id=acc, name=acc, annotations={"organism": species, "molecule_type": "DNA"})
    rec.features.append(SeqFeature(FeatureLocation(1100, 1400, 1), type="CDS"))
    SeqIO.write(rec, path, "genbank")
    return path


def write_submission(tmp_path, entries, organelle="mitochondrion"):
    """entries: [(name, seq, cds_location, organism)]"""
    rows = ["COMMON\tSUBMITTER\t\tab_name\tMishima,T."]
    for name, seq, loc, organism in entries:
        rows.append(f"{name}\tsource\t1..{len(seq)}\torganism\t{organism}")
        rows.append("\t\t\tmol_type\tgenomic DNA")
        if organelle:
            rows.append(f"\t\t\torganelle\t{organelle}")
        rows.append(f"\tCDS\t{loc}\tproduct\tcytochrome c oxidase subunit I")
    ann = tmp_path / "s.ann"
    ann.write_text("\n".join(rows) + "\n")
    fasta = tmp_path / "s.fasta"
    fasta.write_text("".join(f">{n}\n{s}\n//\n" for n, s, _, _ in entries))
    return ann, fasta


def refs_args(tmp_path):
    return ["--refs", str(write_ref(tmp_path / "c.gb", "NC_000001", "Apis cerana", CERANA)),
            str(write_ref(tmp_path / "m.gb", "NC_000002", "Apis mellifera", MELLIFERA))]


@needs_blast
def test_end_to_end_with_refs_tsv_and_rules(tmp_path, capsys):
    q = CERANA[1000:1600]
    ann, fasta = write_submission(tmp_path, [("e1", q, "101..400", "Apis cerana"),
                                             ("e2", q, "104..400", "Apis cerana")])
    rules = tmp_path / "r.json"
    code = cli.main([str(ann), str(fasta), *refs_args(tmp_path), "--tsv", "--rules", str(rules)])
    assert code == 1
    rows = [l.split("\t") for l in capsys.readouterr().out.splitlines()]
    assert rows[0] == ["entry", "check", "result", "detail", "ref", "proposal"]
    assert [(r[0], r[2]) for r in rows[1:]] == [("e2", "BOUNDARY")]
    got = json.loads(rules.read_text())["rules"]
    assert [(r["action"], r["enabled"], r["replacement"]) for r in got] == [("replace_location", False, "101..400")]


@needs_blast
def test_rules_apply_with_fix_edit_after_enabling(tmp_path, capsys):
    q = CERANA[1000:1600]
    ann, fasta = write_submission(tmp_path, [("e2", q, "104..400", "Apis cerana")])
    rules = tmp_path / "r.json"
    cli.main([str(ann), str(fasta), *refs_args(tmp_path), "--rules", str(rules)])
    data = json.loads(rules.read_text())
    for r in data["rules"]:
        r["enabled"] = True
    rules.write_text(json.dumps(data))
    assert edit.main([str(rules), str(ann), "--outdir", str(tmp_path / "out")]) == 0
    assert "\tCDS\t101..400\t" in (tmp_path / "out" / "s.ann").read_text()


@needs_blast
def test_misidentified_entry_is_flagged(tmp_path, capsys):
    ann, fasta = write_submission(tmp_path, [("e1", CERANA[1000:1600], "101..400", "Apis mellifera")])
    assert cli.main([str(ann), str(fasta), *refs_args(tmp_path)]) == 1
    assert "SPECIES_MISMATCH" in capsys.readouterr().out


def test_no_organelle_entries_exits_0(tmp_path, capsys):
    ann, fasta = write_submission(tmp_path, [("e1", CERANA[:300], "1..300", "Apis cerana")], organelle=None)
    assert cli.main([str(ann), str(fasta)]) == 0
    assert "no organelle entries" in capsys.readouterr().out


def test_blastn_missing_exits_2(tmp_path, monkeypatch, capsys):
    ann, fasta = write_submission(tmp_path, [("e1", CERANA[:300], "1..300", "Apis cerana")])
    monkeypatch.setattr(align, "blast_available", lambda: False)
    assert cli.main([str(ann), str(fasta)]) == 2
    assert "blastn" in capsys.readouterr().err


@needs_blast
def test_network_error_exits_2(tmp_path, monkeypatch, capsys):
    ann, fasta = write_submission(tmp_path, [("e1", CERANA[:300], "1..300", "Apis cerana")])

    def boom(self, *a, **k):
        raise refs.ReferenceError("NCBI E-utilities not reachable: down")
    monkeypatch.setattr(refs.NcbiClient, "search", boom)
    assert cli.main([str(ann), str(fasta), "--cache", str(tmp_path / "cache")]) == 2
    assert "not reachable" in capsys.readouterr().err


@needs_blast
def test_send_species_only_when_asked(tmp_path, monkeypatch, capsys):
    ann, fasta = write_submission(tmp_path, [("e1", CERANA[:300], "1..300", "Apis cerana")])
    seen = []

    def spy(self, genus, organelle, species=None):
        seen.append((genus, organelle, species))
        return []
    monkeypatch.setattr(refs.NcbiClient, "search", spy)
    assert cli.main([str(ann), str(fasta), "--cache", str(tmp_path / "c")]) == 0
    assert cli.main([str(ann), str(fasta), "--cache", str(tmp_path / "c"), "--send-species"]) == 0
    assert seen == [("Apis", "mitochondrion", None), ("Apis", "mitochondrion", "Apis cerana")]
    assert "NO_REFERENCE" in capsys.readouterr().out


def test_pick_keeps_every_claimed_species():
    s = [refs.RefSummary(f"NC_{i}", f"Apis s{i}", i, "") for i in range(5)]
    picked, dropped = refs.pick_one_per_species(s, claimed={"Apis s3", "Apis s4"}, max_species=2)
    assert {r.species for r in picked} == {"Apis s3", "Apis s4"} and dropped == 3
