"""nuc.fix.partial_ends on synthetic submissions."""
from ddbj_submission_tools.nuc.ann_fix.annfile import AnnFile
from ddbj_submission_tools.nuc.fix import partial_ends as pe
from ddbj_submission_tools.nuc.mssann import reverse_complement

PRE, POST = "C" * 10, "C" * 20


def cds_seq(start="ATG", stop="TAA", body="GCC" * 8):
    return start + body + stop          # 30 nt: positions 11..40 between PRE and POST


FWD = PRE + cds_seq() + POST            # CDS 11..40, 60 nt
REV = reverse_complement(FWD)           # the same CDS at complement(21..50)


def run(tmp_path, features, seq=FWD, name="e1"):
    rows = [[name, "source", f"1..{len(seq)}", "organism", "Examplea testa"],
            ["", "", "", "mol_type", "genomic DNA"]]
    for key, loc, *quals in features:
        quals = quals or [("product", "hypothetical protein")]
        for i, (k, v) in enumerate(quals):
            rows.append(["", key if i == 0 else "", loc if i == 0 else "", k, v])
    path = tmp_path / "s.ann"
    path.write_text("".join("\t".join(r) + "\n" for r in rows))
    ann = AnnFile.read(path)
    changes = pe.fix(ann, {name: seq})
    locs = [(f.key, f.location) for f in ann.features() if f.key != "source"]
    return locs, changes


# -- CDS -------------------------------------------------------------------

def test_complete_cds_is_unchanged(tmp_path):
    locs, changes = run(tmp_path, [("CDS", "11..40")])
    assert locs == [("CDS", "11..40")] and changes == []


def test_cds_without_start_codon_gets_5prime_partial(tmp_path):
    locs, _ = run(tmp_path, [("CDS", "11..40")], seq=PRE + cds_seq(start="GCC") + POST)
    assert locs == [("CDS", "<11..40")]


def test_cds_without_stop_codon_gets_3prime_partial(tmp_path):
    locs, _ = run(tmp_path, [("CDS", "11..40")], seq=PRE + cds_seq(stop="GCC") + POST)
    assert locs == [("CDS", "11..>40")]


def test_cds_length_not_multiple_of_3_has_no_stop(tmp_path):
    locs, _ = run(tmp_path, [("CDS", "11..41")])
    assert locs == [("CDS", "11..>41")]


def test_transl_table_changes_start_and_stop_codons(tmp_path):
    # vertebrate mitochondrial code (2): ATA starts, AGA stops; neither does in table 1
    seq = PRE + cds_seq(start="ATA", stop="AGA") + POST
    locs, _ = run(tmp_path, [("CDS", "11..40", ("transl_table", "2"))], seq=seq)
    assert locs == [("CDS", "11..40")]
    locs, _ = run(tmp_path, [("CDS", "11..40")], seq=seq)      # no transl_table -> table 1
    assert locs == [("CDS", "<11..>40")]


def test_complement_cds(tmp_path):
    locs, _ = run(tmp_path, [("CDS", "complement(21..50)")], seq=REV)
    assert locs == [("CDS", "complement(21..50)")]
    bad = reverse_complement(PRE + cds_seq(start="GCC") + POST)
    locs, _ = run(tmp_path, [("CDS", "complement(21..50)")], seq=bad)
    assert locs == [("CDS", "complement(21..>50)")]


def test_join_cds(tmp_path):
    seq = PRE + cds_seq(start="GCC")[:10] + "T" * 10 + cds_seq(start="GCC")[10:] + POST
    locs, _ = run(tmp_path, [("CDS", "join(11..20,31..50)")], seq=seq)
    assert locs == [("CDS", "join(<11..20,31..50)")]
    locs, _ = run(tmp_path, [("CDS", "complement(join(21..40,51..60))")], seq=reverse_complement(seq))
    assert locs == [("CDS", "complement(join(21..40,51..>60))")]


def test_codon_start_2_means_5prime_partial(tmp_path):
    locs, _ = run(tmp_path, [("CDS", "10..40", ("codon_start", "2"))], seq="C" * 9 + "G" + cds_seq() + POST)
    assert locs == [("CDS", "<10..40")]


def test_stop_completed_by_transl_except_is_kept(tmp_path):
    seq = PRE + cds_seq()[:-2] + POST              # 11..38 ends with a lone T (TAA by polyadenylation)
    q = [("transl_except", "(pos:38,aa:TERM)"), ("transl_table", "2")]
    locs, _ = run(tmp_path, [("CDS", "11..38", *q)], seq=seq)
    assert locs == [("CDS", "11..38")]


def test_ambiguous_codon_is_reported_not_changed(tmp_path):
    locs, changes = run(tmp_path, [("CDS", "11..40")], seq=PRE + cds_seq(start="NNN") + POST)
    assert locs == [("CDS", "11..40")]
    assert [c.new for c in changes] == [None] and "ambiguous" in changes[0].reason


def test_pseudo_cds_is_skipped(tmp_path):
    locs, _ = run(tmp_path, [("CDS", "11..40", ("pseudo", ""))], seq=PRE + cds_seq(start="GCC") + POST)
    assert locs == [("CDS", "11..40")]


def test_existing_partial_is_kept(tmp_path):
    locs, changes = run(tmp_path, [("CDS", "<11..40")], seq=PRE + cds_seq(start="GCC") + POST)
    assert locs == [("CDS", "<11..40")] and changes == []


# -- mRNA ------------------------------------------------------------------

def test_mrna_ends_shared_with_cds_become_partial(tmp_path):
    locs, _ = run(tmp_path, [("mRNA", "1..40"), ("CDS", "11..40"),
                             ("mRNA", "11..60"), ("mRNA", "1..60")])
    assert locs == [("mRNA", "1..>40"), ("CDS", "11..40"), ("mRNA", "<11..60"), ("mRNA", "1..60")]


def test_mrna_equal_to_cds_is_partial_at_both_ends(tmp_path):
    locs, _ = run(tmp_path, [("mRNA", "11..40"), ("CDS", "11..40")])
    assert locs == [("mRNA", "<11..>40"), ("CDS", "11..40")]


def test_complement_mrna(tmp_path):
    locs, _ = run(tmp_path, [("mRNA", "complement(1..50)"), ("CDS", "complement(21..50)")], seq=REV)
    assert locs == [("mRNA", "complement(1..>50)"), ("CDS", "complement(21..50)")]


def test_joined_mrna(tmp_path):
    seq = PRE + cds_seq()[:10] + "T" * 10 + cds_seq()[10:] + POST
    locs, _ = run(tmp_path, [("mRNA", "join(1..20,31..60)"), ("CDS", "join(11..20,31..50)")], seq=seq)
    assert locs[0] == ("mRNA", "join(1..20,31..60)")
    locs, _ = run(tmp_path, [("mRNA", "join(1..20,31..50)"), ("CDS", "join(11..20,31..50)")], seq=seq)
    assert locs[0] == ("mRNA", "join(1..20,31..>50)")


def test_mrna_follows_cds_made_partial(tmp_path):
    locs, _ = run(tmp_path, [("mRNA", "11..60"), ("CDS", "11..40")], seq=PRE + cds_seq(start="GCC") + POST)
    assert locs == [("mRNA", "<11..60"), ("CDS", "<11..40")]


def test_cds_on_other_strand_does_not_touch_mrna(tmp_path):
    locs, _ = run(tmp_path, [("mRNA", "complement(11..60)"), ("CDS", "11..40")])
    assert locs[0] == ("mRNA", "complement(11..60)")


# -- CLI -------------------------------------------------------------------

def test_cli_outdir_changes_only_the_location(tmp_path):
    ann = tmp_path / "s.ann"
    ann.write_text("e1\tsource\t1..60\torganism\tExamplea testa\n"
                   "\tmRNA\t11..60\tproduct\tx\n\tCDS\t11..40\tproduct\tx\n")
    fasta = tmp_path / "s.fasta"
    fasta.write_text(">e1\n" + PRE + cds_seq(start="GCC") + POST + "\n")
    assert pe.main([str(ann), str(fasta), "--dry-run"]) == 0
    out = tmp_path / "out"
    assert pe.main([str(ann), str(fasta), "--outdir", str(out)]) == 0
    assert (out / "s.ann").read_text() == ("e1\tsource\t1..60\torganism\tExamplea testa\n"
                                           "\tmRNA\t<11..60\tproduct\tx\n\tCDS\t<11..40\tproduct\tx\n")
