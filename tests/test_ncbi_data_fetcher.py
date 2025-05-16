import pytest
import logging
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import Entrez, SeqIO

from deeplotyper.extras.ncbi_data_fetcher import NCBIGeneSequenceFetcher

# Simple context-manager stub to fake Entrez handles
class DummyCM:
    def __init__(self, handle):
        self.handle = handle

    def __enter__(self):
        return self.handle

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

@pytest.fixture(autouse=True)
def clear_caches():
    # Clear lru_caches between tests
    NCBIGeneSequenceFetcher._efetch_records.cache_clear()
    NCBIGeneSequenceFetcher._fetch_linked_record_ids.cache_clear()
    yield

def test_init_default_email_and_warning(caplog):
    caplog.set_level(logging.WARNING)
    fetcher = NCBIGeneSequenceFetcher()
    # Default organism, assembly, and mapping
    assert fetcher.organism == "Homo sapiens"
    assert fetcher.assembly == "GRCh38.p14"
    assert fetcher._assembly_acc == NCBIGeneSequenceFetcher._ASSEMBLY_MAPPING["GRCh38.p14"]
    # Entrez.email was set and warning emitted
    assert Entrez.email == "anonymous@example.com"
    assert "Entrez email not provided" in caplog.text

def test_invalid_assembly_raises():
    with pytest.raises(ValueError) as exc:
        NCBIGeneSequenceFetcher(assembly="INVALID")
    msg = str(exc.value)
    assert "Assembly 'INVALID' is not supported" in msg
    assert "Choose one of" in msg

def test_normalize_genomic_info():
    raw = {
        "AssemblyAccVer": "ASM",
        "AnnotationRelease": "REL",
        "ChrAccVer": "CHR1",
        "ChrStart": "5",
        "ChrStop": "15",
    }
    norm = NCBIGeneSequenceFetcher._normalize_genomic_info(raw)
    assert norm == {
        "assembly": "ASM",
        "annotation_release": "REL",
        "seq_region_name": None,
        "seq_region_accession": "CHR1",
        "start": 5,
        "end": 15,
    }

def test_compute_exon_arrangement_empty():
    arr = NCBIGeneSequenceFetcher._compute_exon_arrangement({})
    assert arr == {"canonical": {}, "alternatives": {}}

def test_compute_exon_arrangement_picks_longest():
    exon_info = {
        "T1": [{"exon_number": 1, "start": 1, "end": 4, "strand": 1, "sequence": "ATGC"}],  # length 4
        "T2": [{"exon_number": 1, "start": 1, "end": 10, "strand": 1, "sequence": "ATGCGGTTAA"}],  # length 10
    }
    arr = NCBIGeneSequenceFetcher._compute_exon_arrangement(exon_info)
    # T2 should be canonical
    assert list(arr["canonical"].keys()) == ["T2"]
    assert arr["canonical"]["T2"] == [1]
    assert arr["alternatives"] == {"T1": [1]}

def test_select_genomic_region_info_matching():
    fetcher = NCBIGeneSequenceFetcher()
    summary = {
        "DocumentSummarySet": {
            "DocumentSummary": [
                {
                    "LocationHist": [
                        {
                            "AssemblyAccVer": fetcher._assembly_acc,
                            "ChrAccVer": "ACC_MATCH",
                            "ChrStart": "100",
                            "ChrStop": "200",
                        }
                    ],
                    "GenomicInfo": [
                        {
                            "ChrAccVer": "ACC_FALLBACK",
                            "ChrStart": "300",
                            "ChrStop": "400",
                        }
                    ],
                }
            ]
        }
    }
    region = fetcher._select_genomic_region_info(summary)
    assert region["ChrAccVer"] == "ACC_MATCH"
    assert region["ChrStart"] == "100"
    assert region["ChrStop"] == "200"

def test_select_genomic_region_info_fallback_and_warning(caplog):
    caplog.set_level(logging.WARNING)
    fetcher = NCBIGeneSequenceFetcher()
    summary = {
        "DocumentSummarySet": {
            "DocumentSummary": [
                {
                    "LocationHist": [
                        {
                            "AssemblyAccVer": "OTHER",
                            "ChrAccVer": "ACC1",
                            "ChrStart": "50",
                            "ChrStop": "60",
                        }
                    ],
                    "GenomicInfo": [
                        {
                            "ChrAccVer": "ACC_FB",
                            "ChrStart": "70",
                            "ChrStop": "80",
                        }
                    ],
                }
            ]
        }
    }
    region = fetcher._select_genomic_region_info(summary)
    assert region["ChrAccVer"] == "ACC_FB"
    assert "Assembly 'GRCh38.p14' not found in LocationHist" in caplog.text

def test_select_genomic_region_info_no_info_raises():
    fetcher = NCBIGeneSequenceFetcher()
    summary = {
        "DocumentSummarySet": {"DocumentSummary": [{"LocationHist": [], "GenomicInfo": []}]}
    }
    with pytest.raises(ValueError):
        fetcher._select_genomic_region_info(summary)

def test_fetch_gene_id_success(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    # esearch context stub
    monkeypatch.setattr(Entrez, "esearch", lambda **kwargs: DummyCM(None))
    # read returns one ID
    monkeypatch.setattr(Entrez, "read", lambda handle: {"IdList": ["12345"]})
    gid = fetcher.fetch_gene_id("MYGENE")
    assert gid == "12345"

def test_fetch_gene_id_not_found(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    monkeypatch.setattr(Entrez, "esearch", lambda **kwargs: DummyCM(None))
    monkeypatch.setattr(Entrez, "read", lambda handle: {"IdList": []})
    with pytest.raises(ValueError) as exc:
        fetcher.fetch_gene_id("MYGENE")
    assert "No gene found for 'MYGENE'" in str(exc.value)

def test_fetch_gene_summary(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    dummy = {"summary": True}
    monkeypatch.setattr(Entrez, "esummary", lambda **kwargs: DummyCM(None))
    monkeypatch.setattr(Entrez, "read", lambda handle: dummy)
    out = fetcher.fetch_gene_summary("123")
    assert out is dummy

def test_fetch_genomic_sequence(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    dummy_rec = SeqRecord(Seq("ATGC"), id="acc1")
    monkeypatch.setattr(Entrez, "efetch", lambda **kwargs: DummyCM(None))
    monkeypatch.setattr(SeqIO, "read", lambda handle, fmt: dummy_rec)
    rec = fetcher.fetch_genomic_sequence("ACC1", 0, 3)
    assert rec is dummy_rec

def test_fetch_cdna_and_protein_sequences(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    dummy_cdna = [SeqRecord(Seq("C"), id="r1")]
    dummy_prot = [SeqRecord(Seq("P"), id="p1")]
    monkeypatch.setattr(fetcher, "_fetch_linked_record_ids", lambda gid, db, linkname: ["r1"] if db == "nuccore" else ["p1"])
    monkeypatch.setattr(fetcher, "_fetch_fasta_records", lambda ids, db: dummy_cdna if db == "nuccore" else dummy_prot)
    assert fetcher.fetch_cdna_sequences("gid") == dummy_cdna
    assert fetcher.fetch_protein_sequences("gid") == dummy_prot

def test_fetch_all_ncbi_sequences_handles_no_mane(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    dummy_g = SeqRecord(Seq("G"), id="g")
    dummy_c = [SeqRecord(Seq("C"), id="c")]
    dummy_p = [SeqRecord(Seq("P"), id="p")]
    dummy_summary = {}
    dummy_region = {"ChrAccVer": "ACC", "ChrStart": "1", "ChrStop": "2"}
    # Patch all steps
    monkeypatch.setattr(fetcher, "fetch_gene_id", lambda name: "gid")
    monkeypatch.setattr(fetcher, "fetch_gene_summary", lambda gid: dummy_summary)
    monkeypatch.setattr(fetcher, "_select_genomic_region_info", lambda sum: dummy_region)
    monkeypatch.setattr(fetcher, "fetch_genomic_sequence", lambda a, s, e: dummy_g)
    monkeypatch.setattr(fetcher, "fetch_cdna_sequences", lambda gid: dummy_c)
    monkeypatch.setattr(fetcher, "fetch_protein_sequences", lambda gid: dummy_p)
    monkeypatch.setattr(fetcher, "get_exon_info_ncbi", lambda gid: {"e": []})
    monkeypatch.setattr(fetcher, "get_exon_arrangement_ncbi", lambda gid: {"a": []})
    # Simulate no MANE
    monkeypatch.setattr(fetcher, "fetch_mane_select_cdna", lambda gid: (_ for _ in ()).throw(ValueError))
    monkeypatch.setattr(fetcher, "fetch_mane_select_protein", lambda gid: (_ for _ in ()).throw(ValueError))

    raw = fetcher.fetch_all_ncbi_sequences("MYGENE")
    assert raw["genomic"] is dummy_g
    assert raw["cdna"] is dummy_c
    assert raw["protein"] is dummy_p
    assert raw["mane_cdna"] is None
    assert raw["mane_protein"] is None
    assert raw["exon_info"] == {"e": []}
    assert raw["region_info"] == dummy_region

def test_fetch_all_transforms_raw_to_simple(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    raw = {
        "genomic": SeqRecord(Seq("XYZ"), id="g"),
        "region_info": {
            "AssemblyAccVer": "A",
            "AnnotationRelease": "R",
            "ChrAccVer": "C",
            "ChrStart": "10",
            "ChrStop": "20",
        },
        "cdna": [SeqRecord(Seq("C1"), id="c1")],
        "protein": [SeqRecord(Seq("P1"), id="p1")],
        "mane_cdna": None,
        "mane_protein": None,
        "exon_info": {"e": []},
        "exon_arrangements": {"a": []},
    }
    monkeypatch.setattr(fetcher, "fetch_all_ncbi_sequences", lambda name: raw)
    out = fetcher.fetch_all("GENE")
    assert out["genomic_sequence"] == "XYZ"
    assert out["genomic_info"] == NCBIGeneSequenceFetcher._normalize_genomic_info(raw["region_info"])
    assert out["cdna"] == {"c1": "C1"}
    assert out["protein"] == {"p1": "P1"}
    assert out["ccds"] == {}
    assert out["mane_cdna"] is None
    assert out["mane_protein"] is None
    assert out["exon_info"] == {"e": []}
    assert out["exon_arrangements"] == {"a": []}

def test_get_exon_info_ncbi(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    # build a SeqRecord with one exon feature
    rec = SeqRecord(Seq("ATGC"), id="T1")
    feat = SeqFeature(FeatureLocation(0, 4, strand=1), type="exon", qualifiers={"number": ["2"]})
    rec.features = [feat]
    # stub out linked IDs and genbank fetch (accept linkname!)
    monkeypatch.setattr(fetcher, "_fetch_linked_record_ids", lambda gid, db, linkname: ["T1"])
    monkeypatch.setattr(fetcher, "_fetch_genbank_records", lambda ids, db: [rec])
    info = fetcher.get_exon_info_ncbi("gid")
    assert "T1" in info
    entry = info["T1"][0]
    assert entry["exon_number"] == 2
    assert entry["start"] == 1
    assert entry["end"] == 4
    assert entry["strand"] == 1
    assert entry["sequence"] == "ATGC"

def test_get_exon_arrangement_ncbi(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    info = {
        "T1": [{"exon_number": 1, "start": 1, "end": 4, "strand": 1, "sequence": "ATGC"}],
        "T2": [{"exon_number": 1, "start": 1, "end": 10, "strand": 1, "sequence": "ATGCGGTTAA"}],
    }
    monkeypatch.setattr(fetcher, "get_exon_info_ncbi", lambda gid: info)
    arr = fetcher.get_exon_arrangement_ncbi("gid")
    assert "canonical" in arr and "alternatives" in arr
    # T2 has the longer total length
    assert list(arr["canonical"].keys()) == ["T2"]
    assert arr["alternatives"] == {"T1": [1]}

def test_get_mane_select_transcript_id_success(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    rec1 = SeqRecord(Seq("A"), id="R1")
    rec1.annotations = {"keywords": ["whatever", "MANE Select"]}
    rec2 = SeqRecord(Seq("B"), id="R2")
    rec2.annotations = {"keywords": []}
    # stub linked IDs (accept linkname!)
    monkeypatch.setattr(fetcher, "_fetch_linked_record_ids", lambda gid, db, linkname: ["R1", "R2"])
    monkeypatch.setattr(fetcher, "_fetch_genbank_records", lambda ids, db: [rec1, rec2])
    tid = fetcher.get_mane_select_transcript_id("gid")
    assert tid == "R1"

def test_get_mane_select_transcript_id_none(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    rec = SeqRecord(Seq("A"), id="R1")
    rec.annotations = {"keywords": []}
    # stub linked IDs (accept linkname!)
    monkeypatch.setattr(fetcher, "_fetch_linked_record_ids", lambda gid, db, linkname: ["R1"])
    monkeypatch.setattr(fetcher, "_fetch_genbank_records", lambda ids, db: [rec])
    with pytest.raises(ValueError):
        fetcher.get_mane_select_transcript_id("gid")

def test_get_mane_select_transcript_id_multiple(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    rec1 = SeqRecord(Seq("A"), id="R1"); rec1.annotations = {"keywords": ["MANE"]}
    rec2 = SeqRecord(Seq("B"), id="R2"); rec2.annotations = {"keywords": ["MANE"]}
    # stub linked IDs (accept linkname!)
    monkeypatch.setattr(fetcher, "_fetch_linked_record_ids", lambda gid, db, linkname: ["R1", "R2"])
    monkeypatch.setattr(fetcher, "_fetch_genbank_records", lambda ids, db: [rec1, rec2])
    with pytest.raises(ValueError) as exc:
        fetcher.get_mane_select_transcript_id("gid")
    assert "Multiple MANE Select transcripts" in str(exc.value)

def test_fetch_mane_select_cdna(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    rec = SeqRecord(Seq("AAA"), id="TID")
    monkeypatch.setattr(fetcher, "get_mane_select_transcript_id", lambda gid: "TID")
    # _efetch_records stub now accepts rettype and parser keywords
    monkeypatch.setattr(
        NCBIGeneSequenceFetcher,
        "_efetch_records",
        lambda self, ids, db, rettype, parser: (rec,)
    )
    out = fetcher.fetch_mane_select_cdna("gid")
    assert out is rec

def test_fetch_mane_select_protein_success(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    cdna = SeqRecord(Seq(""), id="CDNA")
    feat = SeqFeature(FeatureLocation(0, 1), type="CDS", qualifiers={"protein_id": ["PID"]})
    cdna.features = [feat]
    monkeypatch.setattr(fetcher, "fetch_mane_select_cdna", lambda gid: cdna)
    prot = SeqRecord(Seq("PPP"), id="PID")
    monkeypatch.setattr(fetcher, "_fetch_fasta_records", lambda ids, db: [prot])
    out = fetcher.fetch_mane_select_protein("gid")
    assert out is prot

def test_fetch_mane_select_protein_no_cds(monkeypatch):
    fetcher = NCBIGeneSequenceFetcher()
    cdna = SeqRecord(Seq(""), id="CDNA")
    cdna.features = []  # no CDS
    monkeypatch.setattr(fetcher, "fetch_mane_select_cdna", lambda gid: cdna)
    with pytest.raises(ValueError):
        fetcher.fetch_mane_select_protein("gid")

def test__efetch_records_and_caching(monkeypatch):
    # Patch Entrez.efetch and SeqIO.parse to produce predictable output
    seqs = [SeqRecord(Seq("X"), id="1"), SeqRecord(Seq("Y"), id="2")]

    def fake_efetch(db, id, rettype, retmode):
        class CM:
            def __enter__(self_inner):
                return "HANDLE"
            def __exit__(self_inner, a, b, c):
                pass
        return CM()

    def fake_parse(handle, parser):
        assert handle == "HANDLE"
        return seqs

    monkeypatch.setattr(Entrez, "efetch", fake_efetch)
    monkeypatch.setattr(SeqIO, "parse", fake_parse)

    # First call populates the cache
    out1 = NCBIGeneSequenceFetcher._efetch_records(
        ("1", "2"), db="db", rettype="gb", parser="genbank"
    )
    # Second identical call should hit cache (same object)
    out2 = NCBIGeneSequenceFetcher._efetch_records(
        ("1", "2"), db="db", rettype="gb", parser="genbank"
    )
    assert out1 is out2
    assert len(out1) == 2 and out1[0].id == "1"

def test_fetch_fasta_and_genbank_wrappers(monkeypatch):
    seqtuple = (SeqRecord(Seq("A"), id="1"),)
    # Stub out the cached function
    monkeypatch.setattr(
        NCBIGeneSequenceFetcher, "_efetch_records",
        lambda self, ids, db, rettype, parser: seqtuple
    )
    fetcher = NCBIGeneSequenceFetcher()
    fa = fetcher._fetch_fasta_records(["1"], db="db")
    gb = fetcher._fetch_genbank_records(["1"], db="db")
    assert isinstance(fa, list) and fa[0].id == "1"
    assert isinstance(gb, list) and gb[0].id == "1"

def test__fetch_linked_record_ids(monkeypatch):
    # Simulate Entrez.elink → Entrez.read sequence
    linksets = [
        {"LinkSetDb": [
            {"Link": [{"Id": "10"}, {"Id": "20"}]}
        ]}
    ]
    monkeypatch.setattr(Entrez, "elink", lambda **kwargs: DummyCM(None))
    monkeypatch.setattr(Entrez, "read", lambda handle: linksets)
    ids = NCBIGeneSequenceFetcher._fetch_linked_record_ids(
        "gid", db="protein", linkname="ln"
    )
    assert ids == ["10", "20"]
