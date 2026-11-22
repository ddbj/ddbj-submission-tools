"""nuc.organelle_ref.refs: names, organelle filters, one genome per species, NCBI client (no network)."""
import json
import urllib.error

import pytest

from ddbj_submission_tools.nuc.organelle_ref import refs
from ddbj_submission_tools.nuc.organelle_ref.refs import RefSummary


class FakeGet:
    def __init__(self):
        self.urls = []

    def __call__(self, url):
        self.urls.append(url)
        if "esearch" in url:
            return json.dumps({"esearchresult": {"count": "2", "webenv": "W", "querykey": "1"}}).encode()
        if "esummary" in url:
            return json.dumps({"result": {"uids": ["1", "2"],
                                          "1": {"accessionversion": "NC_014295.1", "organism": "Apis cerana",
                                                "taxid": 7461, "title": "Apis cerana mitochondrion, complete genome"},
                                          "2": {"accessionversion": "MN000001.1", "organism": "Apis mellifera",
                                                "taxid": 7460, "title": "Apis mellifera mitochondrion, complete genome"}
                                          }}).encode()
        if "efetch" in url:
            return b"LOCUS       FETCHED\n//\n"
        raise AssertionError(url)


@pytest.fixture
def fake_get():
    return FakeGet()


def test_species_of_unresolved_names():
    assert refs.species_of("Apis cerana japonica") == "Apis cerana"
    for name in ["Apis sp.", "Apis sp. ABC1", "Apis cf. cerana", "Apis aff. cerana", "Apis cerana x Apis mellifera",
                 "uncultured Apis", "Apis"]:
        assert refs.species_of(name) is None, name
    assert refs.genus_of("uncultured bacterium") is None and refs.genus_of("Apis sp.") == "Apis"


def test_organelle_filter():
    assert refs.organelle_filter("mitochondrion") == "mitochondrion"
    assert refs.organelle_filter("mitochondrion:kinetoplast") == "mitochondrion"
    assert refs.organelle_filter("plastid:chloroplast") == "chloroplast"
    assert refs.organelle_filter("plastid:apicoplast") == "plastid"
    with pytest.raises(ValueError):
        refs.organelle_filter("hydrogenosome")


def test_pick_one_per_species_prefers_refseq_keeps_claimed_and_counts_dropped():
    s = [RefSummary("MN1", "Apis a", 1, ""), RefSummary("NC_1", "Apis a", 1, ""),
         RefSummary("MN2", "Apis b", 2, ""), RefSummary("MN3", "Apis c", 3, ""), RefSummary("MN4", "Apis d", 4, "")]
    picked, dropped = refs.pick_one_per_species(s, claimed="Apis d", max_species=2)
    assert [r.accession for r in picked] == ["NC_1", "MN4"]        # RefSeq species first, claimed always in
    assert dropped == 2


def test_search_sends_genus_only_by_default(fake_get):
    client = refs.NcbiClient(get=fake_get, pause=0)
    got = client.search("Apis", "mitochondrion")
    assert all("cerana" not in u for u in fake_get.urls)
    assert "Apis%5BOrganism%5D" in fake_get.urls[0] and "tool=ddbj-submission-tools" in fake_get.urls[0]
    assert "email" not in "".join(fake_get.urls)
    assert {r.accession for r in got} == {"NC_014295", "MN000001"}
    assert {r.species for r in got} == {"Apis cerana", "Apis mellifera"}


def test_http_error_raises_reference_error():
    def boom(url):
        raise urllib.error.URLError("down")
    with pytest.raises(refs.ReferenceError):
        refs.NcbiClient(get=boom, pause=0).search("Apis", "mitochondrion")


def test_bad_json_raises_reference_error():
    with pytest.raises(refs.ReferenceError):
        refs.NcbiClient(get=lambda url: b"<html>busy</html>", pause=0).search("Apis", "mitochondrion")


def test_fetch_uses_cache(tmp_path, fake_get):
    (tmp_path / "NC_014295.gb").write_text("LOCUS cached\n//\n")
    paths = refs.NcbiClient(get=fake_get, pause=0).fetch(["NC_014295", "MN000001"], tmp_path)
    assert [p.name for p in paths] == ["NC_014295.gb", "MN000001.gb"]
    assert sum("efetch" in u for u in fake_get.urls) == 1
    assert (tmp_path / "MN000001.gb").read_text().startswith("LOCUS       FETCHED")


def test_load_refs_reads_species_and_taxid(tmp_path):
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    rec = SeqRecord(Seq("ACGT" * 50), id="NC_000001.1", name="NC_000001",
                    annotations={"organism": "Apis cerana japonica", "molecule_type": "DNA"})
    rec.features.append(SeqFeature(FeatureLocation(0, 200, 1), type="source", qualifiers={"db_xref": ["taxon:7461"]}))
    p = tmp_path / "r.gb"
    SeqIO.write(rec, p, "genbank")
    got = refs.load_refs([p])
    summary, record = got["NC_000001"]
    assert (summary.species, summary.taxid, summary.refseq) == ("Apis cerana", 7461, True)
    assert len(record.seq) == 200


def test_species_of_more_unresolved_forms():
    for name in ["Apis cf cerana", "Apis nr. cerana", "Apis near cerana", "Apis gr. cerana", "Apis spp.", "Apis aff cerana"]:
        assert refs.species_of(name) is None, name


def test_fetch_error_body_is_not_cached(tmp_path):
    client = refs.NcbiClient(get=lambda url: b'{"error":"API rate limit exceeded"}', pause=0)
    with pytest.raises(refs.ReferenceError):
        client.fetch(["NC_999999"], tmp_path)
    assert not (tmp_path / "NC_999999.gb").exists()


def test_search_species(fake_get):
    got = refs.NcbiClient(get=fake_get, pause=0).search_species("Apis cerana", "mitochondrion")
    assert any("Apis+cerana%5BOrganism%5D" in u or "Apis%20cerana%5BOrganism%5D" in u for u in fake_get.urls)
    assert {r.taxid for r in got} == {7461, 7460}
