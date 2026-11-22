"""nuc.check.assembly_gap and nuc.fix.assembly_gap on synthetic submissions."""
from ddbj_submission_tools.nuc.check import assembly_gap as check
from ddbj_submission_tools.nuc.fix import assembly_gap as fix


def seq(*parts):
    """('A', 50) -> 50 A; ('N', 100) -> 100 N."""
    return "".join(base * n for base, n in parts)


def write(tmp_path, ann_rows, seqs):
    ann = tmp_path / "s.ann"
    ann.write_text("".join("\t".join(r) + "\n" for r in ann_rows))
    fasta = tmp_path / "s.fasta"
    fasta.write_text("".join(f">{k}\n{v}\n//\n" for k, v in seqs.items()))
    return ann, fasta


COMMON = [["COMMON", "SUBMITTER", "", "ab_name", "Mishima,T."]]


def entry(name, length):
    return [[name, "source", f"1..{length}", "organism", "Examplea testa"],
            ["", "", "", "mol_type", "genomic DNA"]]


def gap(loc, est="known", gap_type="within scaffold", evidence="paired-ends"):
    rows = [["", "assembly_gap", loc, "estimated_length", est], ["", "", "", "gap_type", gap_type]]
    if evidence:
        rows.append(["", "", "", "linkage_evidence", evidence])
    return rows


# -- distribution and the known / unknown suggestion -----------------------

def test_runs_shorter_than_min_are_not_gaps():
    s = seq(("A", 20), ("N", 9), ("A", 20), ("n", 10), ("A", 5))
    assert check.n_runs(s, 10) == [(50, 59)]


def test_mostly_multiples_of_100_suggest_unknown():
    est, frac = check.suggest_estimated_length([100, 100, 200, 57])
    assert (est, frac) == ("unknown", 0.75)


def test_varied_lengths_suggest_known():
    est, frac = check.suggest_estimated_length([57, 123, 100, 311, 42])
    assert (est, frac) == ("known", 0.2)


def test_distribution_report_lists_lengths_and_counts(tmp_path, capsys):
    s = seq(("A", 10), ("N", 100), ("A", 10), ("N", 100), ("A", 10), ("N", 57), ("A", 10))
    ann, fasta = write(tmp_path, COMMON + entry("c1", len(s)), {"c1": s})
    assert check.main([str(ann), str(fasta)]) == 1          # 3 runs without assembly_gap
    out = capsys.readouterr().out
    assert "100\t2" in out and "57\t1" in out
    assert "unknown" in out


# -- checks on existing gap features ---------------------------------------

def test_existing_gaps_are_checked(tmp_path, capsys):
    s = seq(("A", 10), ("N", 100), ("A", 10), ("N", 50), ("A", 10), ("N", 30), ("A", 10), ("N", 20), ("A", 10))
    # runs: 11..110, 121..170, 181..210, 221..240
    rows = (COMMON + entry("c1", len(s))
            + gap("11..110", est="unknown")                          # ok
            + gap("121..171")                                        # one base past the run
            + gap("181..210", est="25")                              # estimated_length != 30
            + gap("221..240", evidence=None)                         # within scaffold needs linkage_evidence
            )
    ann, fasta = write(tmp_path, rows, {"c1": s})
    found = check.analyze(check.AnnFile.read(ann), check.parse_fasta_sequences(str(fasta)), 10)
    kinds = sorted((p.kind, p.location) for p in found.problems)
    assert kinds == [("estimated_length", "181..210"), ("linkage_evidence", "221..240"),
                     ("location", "121..171")]
    assert found.uncovered == {}                                     # the 121..170 run is overlapped, not missing


def test_bad_vocabulary_is_reported(tmp_path):
    s = seq(("A", 10), ("N", 100), ("A", 10))
    rows = COMMON + entry("c1", len(s)) + gap("11..110", gap_type="inside scaffold", evidence="pair-ends")
    ann, fasta = write(tmp_path, rows, {"c1": s})
    found = check.analyze(check.AnnFile.read(ann), check.parse_fasta_sequences(str(fasta)), 10)
    assert sorted(p.kind for p in found.problems) == ["gap_type", "linkage_evidence"]


def test_sequence_without_entry_is_a_problem(tmp_path):
    s = seq(("A", 10), ("N", 100), ("A", 10))
    ann, fasta = write(tmp_path, COMMON + entry("c1", len(s)), {"c1": s, "c2": s})
    found = check.analyze(check.AnnFile.read(ann), check.parse_fasta_sequences(str(fasta)), 10)
    assert [(p.kind, p.entry) for p in found.problems] == [("entry", "c2")]


def test_check_exits_0_when_all_gaps_are_annotated(tmp_path, capsys):
    s = seq(("A", 10), ("N", 100), ("A", 10))
    ann, fasta = write(tmp_path, COMMON + entry("c1", len(s)) + gap("11..110", est="unknown"), {"c1": s})
    assert check.main([str(ann), str(fasta)]) == 0


# -- adding the missing assembly_gap features ------------------------------

def two_entries(tmp_path):
    s1 = seq(("A", 10), ("N", 100), ("A", 10), ("N", 200), ("A", 10))     # 11..110, 121..320
    s2 = seq(("A", 10), ("N", 100), ("A", 10))                            # 11..110 (already annotated)
    rows = (COMMON + entry("c1", len(s1)) + [["", "CDS", "1..9", "product", "x"]]
            + entry("c2", len(s2)) + gap("11..110", est="unknown"))
    return write(tmp_path, rows, {"c1": s1, "c2": s2})


def test_fix_adds_only_missing_gaps_at_the_end_of_their_entry(tmp_path, capsys):
    ann, fasta = two_entries(tmp_path)
    before = ann.read_text().splitlines()
    assert fix.main([str(ann), str(fasta), "--linkage-evidence", "paired-ends", "--outdir", str(tmp_path / "o")]) == 0
    after = (tmp_path / "o" / "s.ann").read_text().splitlines()
    at = before.index("\tCDS\t1..9\tproduct\tx") + 1          # after c1's last feature, before c2's source
    added = after[at:at + 6]
    assert added == ["\tassembly_gap\t11..110\testimated_length\tunknown",
                     "\t\t\tgap_type\twithin scaffold",
                     "\t\t\tlinkage_evidence\tpaired-ends",
                     "\tassembly_gap\t121..320\testimated_length\tunknown",
                     "\t\t\tgap_type\twithin scaffold",
                     "\t\t\tlinkage_evidence\tpaired-ends"]
    assert after[:at] + after[at + 6:] == before                      # nothing else changed
    assert "estimated_length: unknown" in capsys.readouterr().err   # auto: both runs are multiples of 100


def test_fix_known_writes_known(tmp_path):
    ann, fasta = two_entries(tmp_path)
    fix.main([str(ann), str(fasta), "--linkage-evidence", "paired-ends", "--estimated-length", "known",
              "--outdir", str(tmp_path / "o")])
    assert "\tassembly_gap\t11..110\testimated_length\tknown\n" in (tmp_path / "o" / "s.ann").read_text()


def test_fix_dry_run_writes_nothing(tmp_path):
    ann, fasta = two_entries(tmp_path)
    before = ann.read_text()
    assert fix.main([str(ann), str(fasta), "--linkage-evidence", "paired-ends", "--dry-run"]) == 0
    assert ann.read_text() == before and not (tmp_path / "s.ann.orig").exists()


def test_fix_requires_linkage_evidence_for_within_scaffold(tmp_path, capsys):
    ann, fasta = two_entries(tmp_path)
    assert fix.main([str(ann), str(fasta), "--dry-run"]) == 2
    assert fix.main([str(ann), str(fasta), "--gap-type", "between scaffolds", "--dry-run"]) == 0


def test_fix_rejects_unknown_vocabulary(tmp_path, capsys):
    ann, fasta = two_entries(tmp_path)
    assert fix.main([str(ann), str(fasta), "--linkage-evidence", "pair-ends", "--dry-run"]) == 2


def test_fix_leaves_runs_overlapped_by_a_wrong_gap(tmp_path, capsys):
    s = seq(("A", 10), ("N", 100), ("A", 10))
    ann, fasta = write(tmp_path, COMMON + entry("c1", len(s)) + gap("11..111"), {"c1": s})
    before = ann.read_text()
    assert fix.main([str(ann), str(fasta), "--linkage-evidence", "paired-ends", "--in-place"]) == 1
    assert ann.read_text() == before                                  # nothing to add; the wrong one is reported


def test_no_runs_gives_no_suggestion(tmp_path, capsys):
    ann, fasta = write(tmp_path, COMMON + entry("c1", 30), {"c1": seq(("A", 30))})
    assert check.main([str(ann), str(fasta)]) == 0
    out = capsys.readouterr().out
    assert "no runs of N" in out and "suggestion" not in out
