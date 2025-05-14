import pytest
from types import SimpleNamespace
from deeplotyper.data_models import (
    HaplotypeEvent
)
from deeplotyper.vcf_haploevents import HaplotypeGroups
import pysam

# Helper to generate a fake pysam.VariantFile replacement
# Helper to generate a fake pysam.VariantFile replacement
def make_fake_variantfile(samples, recs):
    class FakeVariantFile:
        def __init__(self, path):
            self.header = SimpleNamespace(samples=samples)
            self._recs = recs

        def fetch(self, contig, start, end):
            return self._recs

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            pass

    return FakeVariantFile

class FakeCall:
    """Fake genotype call object with a .get('GT')."""
    def __init__(self, gt):
        self._gt = gt
    def get(self, key):
        if key == "GT":
            return self._gt
        return None

class FakeRec:
    """Fake VCF record with pos, ref, alts, and sample_calls."""
    def __init__(self, pos, ref, alts, sample_calls):
        self.pos = pos
        self.ref = ref
        self.alts = alts
        self.samples = sample_calls

def test_build_sequence_various_events():
    ref = "ATCG"
    ev_snv = HaplotypeEvent(pos0=1, ref_allele="T", alt_seq="G")
    ev_ins = HaplotypeEvent(pos0=2, ref_allele="C", alt_seq="GGG")
    ev_del = HaplotypeEvent(pos0=1, ref_allele="TC", alt_seq="")
    hg = HaplotypeGroups({}, ref_seq=ref, start=0)

    assert hg.build_sequence((ev_snv,)) == "AGCG"
    assert hg.build_sequence((ev_ins,)) == "ATGGGG"
    assert hg.build_sequence((ev_del,)) == "AG"
    # out-of-order should sort by pos0
    assert hg.build_sequence((ev_ins, ev_snv)) == "AGGGGG"

def test_build_sequence_clamping():
    ref = "ATCG"
    ev_left = HaplotypeEvent(pos0=-1, ref_allele="A", alt_seq="G")
    ev_right = HaplotypeEvent(pos0=10, ref_allele="G", alt_seq="T")
    hg = HaplotypeGroups({}, ref_seq=ref, start=0)

    assert hg.build_sequence((ev_left,)) == "GATCG"
    assert hg.build_sequence((ev_right,)) == "ATCGT"

def test_materialize_mapping():
    ev1 = HaplotypeEvent(0, "A", "G")
    ev2 = HaplotypeEvent(1, "C", "T")
    mapping = {
        (ev1,): ("S1",),
        (ev2,): ("S2", "S3")
    }
    hg = HaplotypeGroups(mapping, ref_seq="AC", start=0)
    mat = hg.materialize()

    assert mat["GC"] == ("S1",)
    assert mat["AT"] == ("S2", "S3")

def test_from_vcf_missing_samples(monkeypatch):
    FakeVF = make_fake_variantfile(["S1","S2"], [])
    monkeypatch.setattr(pysam, "VariantFile", FakeVF)
    with pytest.raises(ValueError):
        HaplotypeGroups.from_vcf(
            vcf_path="x.vcf", ref_seq="A", contig="chr", start=0,
            samples=["S3"]
        )

def test_from_vcf_default_and_valid(monkeypatch):
    rec = FakeRec(
        pos=1, ref="C", alts=("T",),
        sample_calls={"S1": FakeCall((1,0))}
    )
    FakeVF = make_fake_variantfile(["S1"], [rec])
    monkeypatch.setattr(pysam, "VariantFile", FakeVF)

    hg = HaplotypeGroups.from_vcf(
        vcf_path="x.vcf", ref_seq="ACG", contig="chr", start=0
    )
    mat = hg.materialize()
    assert mat["ATG"] == ("S1",)
    assert mat["ACG"] == ("S1",)

def test_from_vcf_skip_out_of_window(monkeypatch):
    rec = FakeRec(pos=3, ref="C", alts=("T",), sample_calls={"S1": FakeCall((1,))})
    FakeVF = make_fake_variantfile(["S1"], [rec])
    monkeypatch.setattr(pysam, "VariantFile", FakeVF)

    hg = HaplotypeGroups.from_vcf(
        vcf_path="x.vcf", ref_seq="AC", contig="chr", start=0
    )
    assert hg.materialize() == {"AC": ("S1",)}

def test_from_vcf_skip_missing_gt(monkeypatch):
    rec = FakeRec(pos=0, ref="A", alts=("C",), sample_calls={"S1": FakeCall((None,))})
    FakeVF = make_fake_variantfile(["S1"], [rec])
    monkeypatch.setattr(pysam, "VariantFile", FakeVF)

    hg = HaplotypeGroups.from_vcf(
        vcf_path="x.vcf", ref_seq="A", contig="chr", start=0
    )
    assert hg.materialize() == {"A": ("S1",)}

def test_from_vcf_skip_symbolic_and_oob(monkeypatch):
    rec1 = FakeRec(pos=0, ref="A", alts=("<DEL>",), sample_calls={"S1": FakeCall((1,))})
    rec2 = FakeRec(pos=0, ref="A", alts=("G",), sample_calls={"S1": FakeCall((2,))})
    FakeVF = make_fake_variantfile(["S1"], [rec1, rec2])
    monkeypatch.setattr(pysam, "VariantFile", FakeVF)

    hg = HaplotypeGroups.from_vcf(
        vcf_path="x.vcf", ref_seq="A", contig="chr", start=0
    )
    assert hg.materialize() == {"A": ("S1",)}

def test_build_sequence_long_ref_large_indels():
    ref = "ATCGATCGATCG"
    hg = HaplotypeGroups({}, ref_seq=ref, start=0)

    # insertion at beginning
    ev_ins_begin = HaplotypeEvent(0, "A", "GGGG")
    assert hg.build_sequence((ev_ins_begin,)) == "GGGGTCGATCGATCG"

    # deletion at beginning
    ev_del_begin = HaplotypeEvent(0, "ATC", "")
    assert hg.build_sequence((ev_del_begin,)) == "GATCGATCG"

    # insertion in middle
    ev_ins_mid = HaplotypeEvent(5, "T", "CCCCCC")
    assert hg.build_sequence((ev_ins_mid,)) == "ATCGACCCCCCCGATCG"

    # deletion in middle
    ev_del_mid = HaplotypeEvent(2, "CGAT", "")
    assert hg.build_sequence((ev_del_mid,)) == "ATCGATCG"

    # insertion at end
    ev_ins_end = HaplotypeEvent(len(ref)-1, "G", "TTTTT")
    assert hg.build_sequence((ev_ins_end,)) == "ATCGATCGATCTTTTT"

    # deletion at end
    ev_del_end = HaplotypeEvent(len(ref)-4, "ATCG", "")
    assert hg.build_sequence((ev_del_end,)) == "ATCGATCG"

    # combined: begin-insert, mid-delete, end-insert
    ev_del_mid2 = HaplotypeEvent(5, "T", "")
    combined = (ev_ins_begin, ev_del_mid2, ev_ins_end)
    # should yield 'GGGG' + 'TCGA' + '' + 'CGATC' + 'TTTTT'
    assert hg.build_sequence(combined) == "GGGGTCGACGATCTTTTT"