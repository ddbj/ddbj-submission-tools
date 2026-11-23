"""nuc.fix.derivable_features on synthetic submissions."""
from ddbj_submission_tools.nuc.ann_fix.annfile import AnnFile
from ddbj_submission_tools.nuc.fix import derivable_features as df

GENE = [("mRNA", "join(100..200,300..400)"), ("CDS", "join(150..200,300..350)"),
        ("exon", "100..200"), ("intron", "201..299"), ("exon", "300..400"),
        ("5'UTR", "100..149"), ("3'UTR", "351..400")]


def run(tmp_path, features, source=(("organism", "Examplea testa"), ("mol_type", "genomic DNA")),
        common=(), keys=None):
    rows = [["COMMON", key if i == 0 else "", "", k, v] for key, quals in common for i, (k, v) in enumerate(quals)]
    for i, (k, v) in enumerate(source):
        rows.append(["e1" if i == 0 else "", "source" if i == 0 else "", "1..1000" if i == 0 else "", k, v])
    for key, loc, *quals in features:
        quals = quals or [("locus_tag", "EXA_0001")]
        for i, (k, v) in enumerate(quals):
            rows.append(["", key if i == 0 else "", loc if i == 0 else "", k, v])
    path = tmp_path / "s.ann"
    path.write_text("".join("\t".join(r) + "\n" for r in rows))
    ann = AnnFile.read(path)
    results = df.drop(ann, keys or df.DEFAULT_KEYS)
    left = [(f.key, f.location) for f in ann.features() if f.key not in ("source", "DIVISION")]
    return left, results


def test_parts_inside_mrna_with_cds_are_dropped(tmp_path):
    left, results = run(tmp_path, GENE)
    assert left == [("mRNA", "join(100..200,300..400)"), ("CDS", "join(150..200,300..350)")]
    assert sum(r.dropped for r in results) == 5


def test_keys_option_limits_what_is_dropped(tmp_path):
    left, _ = run(tmp_path, GENE, keys={"intron"})
    assert ("intron", "201..299") not in left and ("exon", "100..200") in left


def test_exon_outside_any_mrna_is_kept(tmp_path):
    # e.g. exons of a ncRNA gene: no mRNA to derive them from
    left, results = run(tmp_path, [("ncRNA", "join(100..200,300..400)"), ("exon", "100..200")])
    assert ("exon", "100..200") in left
    assert not results[0].dropped and "not within an mRNA" in results[0].reason


def test_exon_in_an_mrna_intron_is_kept(tmp_path):
    # e.g. an exon of another gene nested in the intron
    left, _ = run(tmp_path, GENE[:2] + [("exon", "220..280"), ("intron", "150..250")])
    assert ("exon", "220..280") in left and ("intron", "150..250") in left


def test_utr_needs_a_cds_in_the_mrna(tmp_path):
    left, _ = run(tmp_path, [("mRNA", "100..400"), ("5'UTR", "100..149")])
    assert ("5'UTR", "100..149") in left


def test_other_strand_mrna_does_not_count(tmp_path):
    left, _ = run(tmp_path, [("mRNA", "complement(100..400)"), ("CDS", "complement(150..350)"),
                             ("exon", "100..400")])
    assert ("exon", "100..400") in left


def test_complement_parts_are_dropped(tmp_path):
    left, _ = run(tmp_path, [("mRNA", "complement(100..400)"), ("CDS", "complement(150..350)"),
                             ("exon", "complement(100..400)"), ("3'UTR", "complement(100..149)")])
    assert left == [("mRNA", "complement(100..400)"), ("CDS", "complement(150..350)")]


def test_feature_with_a_note_is_kept(tmp_path):
    left, results = run(tmp_path, GENE[:2] + [("exon", "100..200", ("locus_tag", "EXA_0001"),
                                               ("note", "alternative exon"))])
    assert ("exon", "100..200") in left and "note" in results[0].reason


def test_allow_qualifier_drops_features_with_it(tmp_path):
    rows = GENE[:2] + [("intron", "201..299", ("note", "source:AUGUSTUS"))]
    path = tmp_path / "a.ann"
    path.write_text("e1\tsource\t1..1000\torganism\tExamplea testa\n"
                    + "".join(f"\t{k}\t{l}\t{q[0][0] if q else 'locus_tag'}\t{q[0][1] if q else 'x'}\n"
                              for k, l, *q in rows))
    ann = AnnFile.read(path)
    assert df.drop(ann, df.DEFAULT_KEYS, ["note"])[0].dropped


def test_organelle_entry_is_skipped(tmp_path):
    left, results = run(tmp_path, GENE, source=(("organism", "Examplea testa"), ("organelle", "mitochondrion")))
    assert len(left) == len(GENE) and results == []


def test_virus_is_skipped(tmp_path):
    left, _ = run(tmp_path, GENE, source=(("organism", "Example virus 1"),))
    assert len(left) == len(GENE)
    left, _ = run(tmp_path, GENE, common=[("DIVISION", [("division", "VRL")])])
    assert len(left) == len(GENE)


def test_organelle_in_common_source_is_skipped(tmp_path):
    left, _ = run(tmp_path, GENE, common=[("source", [("organelle", "plastid:chloroplast")])])
    assert len(left) == len(GENE)


def test_cli_dry_run_and_outdir(tmp_path):
    ann = tmp_path / "s.ann"
    ann.write_text("e1\tsource\t1..1000\torganism\tExamplea testa\n"
                   "\tmRNA\t100..400\tlocus_tag\tEXA_0001\n"
                   "\tCDS\t150..350\tlocus_tag\tEXA_0001\n"
                   "\texon\t100..400\tlocus_tag\tEXA_0001\n"
                   "\t\t\tnumber\t1\n")
    assert df.main([str(ann), "--dry-run"]) == 0
    assert "exon" in ann.read_text()
    out = tmp_path / "out"
    assert df.main([str(ann), "--outdir", str(out)]) == 0
    assert (out / "s.ann").read_text() == ("e1\tsource\t1..1000\torganism\tExamplea testa\n"
                                           "\tmRNA\t100..400\tlocus_tag\tEXA_0001\n"
                                           "\tCDS\t150..350\tlocus_tag\tEXA_0001\n")


def test_all_drops_every_feature_of_the_keys(tmp_path):
    feats = [("misc_RNA", "join(100..200,300..400)"), ("intron", "201..299", ("note", "source:AUGUSTUS")),
             ("mRNA", "join(100..200,300..400)"), ("3'UTR", "180..400")]
    path = tmp_path / "s.ann"
    path.write_text("e1\tsource\t1..1000\torganism\tExamplea testa\n"
                    + "".join(f"\t{k}\t{l}\t{q[0][0] if q else 'locus_tag'}\t{q[0][1] if q else 'x'}\n"
                              for k, l, *q in feats))
    ann = AnnFile.read(path)
    results = df.drop(ann, df.DEFAULT_KEYS, everything=True)
    assert [f.key for f in ann.features()] == ["source", "misc_RNA", "mRNA"]
    assert all(r.dropped for r in results)


def test_all_still_skips_organelles(tmp_path):
    left, _ = run(tmp_path, GENE, source=(("organism", "Examplea testa"), ("organelle", "mitochondrion")))
    assert len(left) == len(GENE)
    path = tmp_path / "s.ann"
    assert df.main([str(path), "--all", "--dry-run"]) == 0
