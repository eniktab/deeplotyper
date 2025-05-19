import pytest
import logging
from typing import Any, Dict, List, Optional, Tuple
from requests.exceptions import HTTPError

from deeplotyper.extras.ensembl_sequence_fetcher import EnsemblGeneDataFetcher


class DummyResponse:
    """A fake requests.Response for testing _get_json."""

    def __init__(self, json_data: Any = None, ok: bool = True):
        self._json = json_data if json_data is not None else {}
        self.ok = ok
        self.raise_called = False

    def raise_for_status(self):
        self.raise_called = True
        if not self.ok:
            raise HTTPError("HTTP error")

    def json(self):
        return self._json


class DummySession:
    """A fake Session whose .get returns a preconfigured DummyResponse."""

    def __init__(self, response: DummyResponse):
        self.response = response
        self.headers: Dict[str, str] = {}
        self.requests: List[Tuple[str, Optional[Dict[str, Any]]]] = []

    def get(self, url: str, params: Optional[Dict[str, Any]] = None):
        self.requests.append((url, params))
        return self.response


@pytest.fixture
def fetcher():
    return EnsemblGeneDataFetcher("MYGENE", species="test_species")


def test_init_sets_fields_and_headers(fetcher, monkeypatch):
    # On init, .gene_symbol, .species and session headers should be set
    f = EnsemblGeneDataFetcher("ABC", species="mus_musculus")
    assert f.gene_symbol == "ABC"
    assert f.species == "mus_musculus"
    # The session is a requests.Session with JSON headers
    assert "Content-Type" in f._session.headers
    assert f._session.headers["Content-Type"] == "application/json"


def test_get_json_success(monkeypatch, fetcher):
    dummy = DummyResponse(json_data={"foo": "bar"}, ok=True)
    # replace the session with our dummy
    monkeypatch.setattr(fetcher, "_session", DummySession(dummy))
    out = fetcher._get_json("/some/endpoint", params={"a": 1})
    # should return the JSON
    assert out == {"foo": "bar"}
    # raise_for_status was called
    assert dummy.raise_called


def test_get_json_http_error(monkeypatch, fetcher):
    dummy = DummyResponse(json_data=None, ok=False)
    monkeypatch.setattr(fetcher, "_session", DummySession(dummy))
    with pytest.raises(HTTPError):
        fetcher._get_json("/bad/endpoint")


def test_fetch_sequence_present(monkeypatch, fetcher):
    # stub _get_json to return a dict with 'seq'
    monkeypatch.setattr(
        fetcher,
        "_get_json",
        lambda p,
        params=None: {
            "seq": "ATGC"})
    assert fetcher._fetch_sequence("/whatever") == "ATGC"


def test_fetch_sequence_missing(monkeypatch, fetcher):
    monkeypatch.setattr(fetcher, "_get_json", lambda p, params=None: {})
    assert fetcher._fetch_sequence("/whatever") == ""


def test_normalize_genomic_info():
    raw = {
        "assembly_name": "ASM01",
        "seq_region_name": "chrX",
        "start": 10,
        "end": 20
    }
    norm = EnsemblGeneDataFetcher._normalize_genomic_info(raw)
    assert norm == {
        "assembly": "ASM01",
        "annotation_release": None,
        "seq_region_name": "chrX",
        "seq_region_accession": None,
        "start": 10,
        "end": 20
    }


def test_get_genomic_sequence(monkeypatch, fetcher):
    lookup = {"seq_region_name": "5", "start": 100, "end": 200, "strand": -1}
    # capture the calls
    calls = []

    def fake_get_json(endpoint, params=None):
        calls.append((endpoint, params))
        return lookup
    monkeypatch.setattr(fetcher, "_get_json", fake_get_json)
    monkeypatch.setattr(fetcher, "_fetch_sequence", lambda path: "GENSEQ")
    out_lookup, out_seq = fetcher.get_genomic_sequence()
    # verify we got back the lookup JSON and the stubbed sequence
    assert out_lookup is lookup
    assert out_seq == "GENSEQ"
    # verify the lookup endpoint was correct
    expected_lookup_ep = f"/lookup/symbol/{fetcher.species}/{fetcher.gene_symbol}"
    assert calls[0] == (expected_lookup_ep, None)


def test_get_ccds_records_all_paths(monkeypatch, fetcher, caplog):
    # 1) Non-CCDS xref is skipped, CCS fetch succeeds, N stripped, HTTPError
    # on xref returns {}
    caplog.set_level(logging.DEBUG)
    xrefs = [
        {"id": "NOT_CCDS1"},
        {"id": "CCDS123"},
        {"id": "CCDS_BAD"}
    ]
    # stub first _get_json for xrefs

    def fake_get_json(endpoint, params=None):
        if endpoint.startswith("/xrefs"):
            return xrefs
        raise AssertionError("unexpected endpoint")
    monkeypatch.setattr(fetcher, "_get_json", fake_get_json)
    # stub _fetch_sequence to raise HTTPError for CCDS_BAD, return seq for
    # CCDS123

    def fake_fetch_sequence(path, params=None):
        if "CCDS_BAD" in path:
            raise HTTPError("fail")
        return "AAAANN"
    monkeypatch.setattr(fetcher, "_fetch_sequence", fake_fetch_sequence)
    records = fetcher.get_ccds_records()
    # should only have CCDS123 and N's stripped
    assert records == {"CCDS123": "AAAA"}
    # verify debug log for skipping non-CCDS
    assert "Skipping non-CCDS xref" in caplog.text

    # 2) HTTPError on xrefs call returns {}
    def bad_json(endpoint, params=None):
        raise HTTPError("no xrefs")
    monkeypatch.setattr(fetcher, "_get_json", bad_json)
    assert fetcher.get_ccds_records() == {}


def test_get_transcript_info_errors_and_success(monkeypatch, fetcher):
    # a) ValueError when no transcripts
    monkeypatch.setattr(
        fetcher,
        "_get_json",
        lambda ep,
        params=None: {
            "Transcript": []})
    with pytest.raises(ValueError):
        fetcher.get_transcript_info()

    # b) Normal case with an explicit canonical
    transcripts = [
        {"id": "TX1", "biotype": "b1", "Translation": {"id": "P1"}, "is_canonical": 0},
        {"id": "TX2", "biotype": "b2", "Translation": {"id": "P2"}, "is_canonical": 1},
    ]
    # first call to _get_json returns the gene_data

    def fake_get_json(endpoint, params=None):
        if "lookup/symbol" in endpoint:
            return {"Transcript": transcripts}
        raise AssertionError("unhandled endpoint")
    monkeypatch.setattr(fetcher, "_get_json", fake_get_json)
    # stub seq fetcher: return cdna or protein based on params

    def fake_fetch_sequence(path, params=None):
        return "PROT" if params and params.get("type") == "protein" else "CDNA"
    monkeypatch.setattr(fetcher, "_fetch_sequence", fake_fetch_sequence)

    info = fetcher.get_transcript_info()
    # should pick TX2 as canonical
    assert "TX2" in info["canonical"]
    assert info["canonical"]["TX2"]["cdna"] == "CDNA"
    # alternative TX1
    assert "TX1" in info["alternatives"]

    # c) fallback when no is_canonical: longest cdna
    transcripts2 = [
        {"id": "A", "biotype": "b", "Translation": {}, "is_canonical": 0},
        {"id": "B", "biotype": "b", "Translation": {}, "is_canonical": 0},
    ]

    def fake_get_json2(endpoint, params=None):
        return {"Transcript": transcripts2}
    monkeypatch.setattr(fetcher, "_get_json", fake_get_json2)
    # stub sequence lengths: A->"AA", B->"BBBB"

    def fake_fetch2(path, params=None):
        return "BBBB" if "B" in path else "AA"
    monkeypatch.setattr(fetcher, "_fetch_sequence", fake_fetch2)

    info2 = fetcher.get_transcript_info()
    # should pick B as canonical
    assert list(info2["canonical"].keys()) == ["B"]
    assert list(info2["alternatives"].keys()) == ["A"]


def test_combine_transcript_info():
    ti = {
        "canonical": {"C1": {"cdna": "AAA", "protein": "PPP"}},
        "alternatives": {"A1": {"cdna": "CCC", "protein": None}}
    }
    cdna_map, prot_map = EnsemblGeneDataFetcher.combine_transcript_info(ti)
    assert cdna_map == {"C1": "AAA", "A1": "CCC"}
    assert prot_map == {"C1": "PPP", "A1": None}


def test_get_exon_info_and_sort(monkeypatch, fetcher):
    # stub get_transcript_info to return a single canonical
    monkeypatch.setattr(
        fetcher,
        "get_transcript_info",
        lambda: {
            "canonical": {
                "T1": {}},
            "alternatives": {}})
    # stub _get_json to return two exons out of order
    exons = [
        {"id": "E1", "seq_region_name": "1", "start": 10, "end": 20, "strand": 1},
        {"id": "E2", "seq_region_name": "1", "start": 5, "end": 9, "strand": 1},
    ]
    monkeypatch.setattr(
        fetcher,
        "_get_json",
        lambda ep,
        params=None: {
            "Exon": exons})
    # stub _fetch_sequence to echo the region string
    monkeypatch.setattr(fetcher, "_fetch_sequence",
                        lambda path: path.split("/")[-1])

    exon_info = fetcher.get_exon_info()
    # entries should be sorted by start: E2 then E1
    entries = exon_info["T1"]
    assert entries[0]["exon_id"] == "E2"
    assert entries[0]["exon_number"] == 1
    # the returned sequence should match the region part
    assert entries[0]["sequence"].startswith("1:5..9:1")


def test_get_exon_arrangement_ensembl(monkeypatch, fetcher):
    # stub both dependencies
    ti = {"canonical": {"C": {}}, "alternatives": {"A": {}}}
    ei = {"C": [{"exon_number": 1}, {"exon_number": 2}],
          "A": [{"exon_number": 3}]}
    monkeypatch.setattr(fetcher, "get_transcript_info", lambda: ti)
    monkeypatch.setattr(fetcher, "get_exon_info", lambda ids=None: ei)

    arr = fetcher.get_exon_arrangement_ensembl()
    assert arr == {"canonical": {"C": [1, 2]}, "alternatives": {"A": [3]}}


def test_fetch_all_ensembl_sequences_and_fetch_all(monkeypatch, fetcher):
    # prepare dummy return for fetch_all_ensembl_sequences
    raw = {
        "gene_record": "GENOME_FAKE",
        "genomic_info": {
            "assembly_name": "X",
            "seq_region_name": "Y",
            "start": 1,
            "end": 2},
        "ccds": {
            "C1": "CCDSSEQ"},
        "cdna": {
            "T1": "CDNASEQ"},
        "protein": {
            "T1": "PROTSEQ"},
        "exon_arrangements": {
            "canonical": {
                "T1": [1]},
            "alternatives": {}},
        "exon_info_canonical": {
            "T1": [
                {
                    "exon_id": "E",
                    "exon_number": 1,
                    "seq_region": "Y",
                    "start": 1,
                    "end": 2,
                    "strand": 1,
                    "sequence": "ESEQ"}]},
    }
    monkeypatch.setattr(fetcher, "fetch_all_ensembl_sequences", lambda: raw)
    # test fetch_all
    out = fetcher.fetch_all()
    # check that fetch_all normalized genomic_info and copied everything else
    assert out["genomic_sequence"] == "GENOME_FAKE"
    assert out["genomic_info"] == {
        "assembly": "X",
        "annotation_release": None,
        "seq_region_name": "Y",
        "seq_region_accession": None,
        "start": 1,
        "end": 2
    }
    assert out["ccds"] == raw["ccds"]
    assert out["cdna"] == raw["cdna"]
    assert out["protein"] == raw["protein"]
    # Ensembl has no MANE, should be None
    assert out["mane_cdna"] is None and out["mane_protein"] is None
    assert out["exon_arrangements"] == raw["exon_arrangements"]
    assert out["exon_info"] == raw["exon_info_canonical"]
