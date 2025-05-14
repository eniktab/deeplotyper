import pytest
from Bio.Seq import Seq

from deeplotyper import HaplotypeRemapper
from deeplotyper.data_models import HaplotypeEvent, NewTranscriptSequences


class DummyResult:
    """Minimal stand-in for TranscriptMappingResult with only the attributes that
    HaplotypeRemapper actually uses."""
    def __init__(self, offset: int, cdna_to_dna_map: dict[int, int], seq_region: str):
        self.offset = offset
        self.cdna_to_dna_map = cdna_to_dna_map
        self.seq_region = seq_region


@pytest.fixture
def simple_mapper():
    # A 6-bp reference with a single transcript that covers the entire slice.
    ref = "ATGAAA"
    mapping = {i + 1: i for i in range(len(ref))}  # cdna positions 1–6 → genome positions 0–5
    res = DummyResult(offset=0, cdna_to_dna_map=mapping, seq_region="chr1")
    return HaplotypeRemapper(reference_genome=ref, transcript_results={"tx1": res})


def test_initial_offset_empty_results():
    remapper = HaplotypeRemapper(reference_genome="XXX", transcript_results={})
    # With no transcripts, offset should default to 0
    assert remapper._offset == 0


def test_empty_haplotype_map(simple_mapper):
    result = simple_mapper.apply_haplotypes({})
    # No haplotypes → empty result
    assert result == {}


def test_no_events(simple_mapper):
    hap_map = {tuple(): tuple()}
    result = simple_mapper.apply_haplotypes(hap_map)
    key = frozenset()
    assert list(result.keys()) == [key]

    out = result[key]
    # 1) genome unchanged
    assert out["mutated_genome"] == "ATGAAA"
    # 2) protein is the translation of the full ORF (no STOP codons → whole seq)
    assert out["genome_protein"] == str(Seq("ATGAAA").translate())
    # 3) identity mapping
    assert out["orig_to_mutated_index"] == {i: i for i in range(6)}
    assert out["mutated_index_to_orig"] == list(range(6))
    # 4) transcript rebuilt unchanged
    tx = out["transcripts"]["tx1"]
    assert isinstance(tx, NewTranscriptSequences)
    assert tx.cdna_sequence == "ATGAAA"
    # 5) codon_map has exactly two codons: ATG→M, AAA→K
    codon_map = tx.codon_map
    assert set(codon_map) == {1, 2}
    assert codon_map[1].codon == "ATG" and codon_map[1].amino_acid == "M"
    assert codon_map[2].codon == "AAA" and codon_map[2].amino_acid == "K"


def test_substitution(simple_mapper):
    # Change G→T at zero-based position 2 (third base)
    ev = HaplotypeEvent(pos0=2, ref_allele="G", alt_seq="T")
    hap_map = {(ev,): tuple()}
    result = simple_mapper.apply_haplotypes(hap_map)
    out = result[frozenset((ev,))]

    assert out["mutated_genome"] == "ATTAAA"
    assert out["genome_protein"] == str(Seq("ATTAAA").translate())
    # mapping unchanged
    assert out["orig_to_mutated_index"][2] == 2
    assert out["mutated_index_to_orig"][2] == 2

    tx = out["transcripts"]["tx1"]
    assert tx.cdna_sequence == "ATTAAA"
    # ATT→I, AAA→K
    assert tx.codon_map[1].amino_acid == "I"
    assert tx.codon_map[2].amino_acid == "K"


def test_deletion(simple_mapper):
    # Delete the first codon "ATG"
    ev = HaplotypeEvent(pos0=0, ref_allele="ATG", alt_seq="")
    hap_map = {(ev,): tuple()}
    result = simple_mapper.apply_haplotypes(hap_map)
    out = result[frozenset((ev,))]

    assert out["mutated_genome"] == "AAA"
    assert out["genome_protein"] == str(Seq("AAA").translate())
    # Positions 0,1,2 are deleted
    for i in (0, 1, 2):
        assert i not in out["orig_to_mutated_index"]
    # The remainder shift down:
    assert out["orig_to_mutated_index"] == {3: 0, 4: 1, 5: 2}
    assert out["mutated_index_to_orig"] == [3, 4, 5]


def test_insertion_length_gt1(simple_mapper):
    # Insert "CCC" (3 bp, to keep codons in‐frame) at pos0=1
    ev = HaplotypeEvent(pos0=1, ref_allele="", alt_seq="CCC")
    hap_map = {(ev,): tuple()}
    result = simple_mapper.apply_haplotypes(hap_map)
    out = result[frozenset((ev,))]

    # Insertion index = pos0 + shift = 1
    assert out["mutated_genome"] == "ACCCTGAAA"
    # Validate protein translation on the 9-bp result
    assert out["genome_protein"] == str(Seq("ACCCTGAAA").translate())

    # orig→mut: base 0 stays at 0; base 1 moves to 1+3=4
    assert out["orig_to_mutated_index"][0] == 0
    assert out["orig_to_mutated_index"][1] == 4

    # The three inserted positions (1,2,3) have no original source
    mt2o = out["mutated_index_to_orig"]
    for idx in (1, 2, 3):
        assert mt2o[idx] is None


def test_multiple_events_order(simple_mapper):
    # Combine a substitution at 0 and an insertion at 2, but give them reversed
    ev_sub = HaplotypeEvent(pos0=0, ref_allele="A", alt_seq="G")
    ev_ins = HaplotypeEvent(pos0=2, ref_allele="", alt_seq="TTT")
    hap_map = {(ev_ins, ev_sub): tuple()}
    result = simple_mapper.apply_haplotypes(hap_map)
    out = result[frozenset((ev_ins, ev_sub))]
    # Should apply sub first (pos0=0), then insertion (pos0=2)
    # 1) ATGAAA → GTGAAA
    # 2) insert "TTT" at idx=2 → GTTTTGAAA
    assert out["mutated_genome"] == "GTTTTGAAA"


def test_wrong_ref_allele(simple_mapper):
    # If the reference allele doesn't match, we should get an AssertionError
    ev_bad = HaplotypeEvent(pos0=1, ref_allele="X", alt_seq="A")
    with pytest.raises(AssertionError):
        simple_mapper.apply_haplotypes({(ev_bad,): tuple()})
