import pytest
from Bio.Align import PairwiseAligner

from deeplotyper import make_aligner, apply_alignment_gaps
from deeplotyper import SequenceCoordinateMapper
from deeplotyper.data_models import (
    RawBase,
    BaseCoordinateMapping,
    TranscriptMappingResult,
)


# ─── Tests for alignment_utils.py ────────────────────────────────────────────────

def test_make_aligner_properties():
    aln = make_aligner(
        mode="global",
        match_score=1.0,
        mismatch_score=-2.0,
        open_gap_score=-5.0,
        extend_gap_score=-2.0,
    )
    assert isinstance(aln, PairwiseAligner)
    assert aln.mode == "global"
    assert aln.match_score == 1.0
    assert aln.mismatch_score == -2.0
    assert aln.open_gap_score == -5.0
    assert aln.extend_gap_score == -2.0


def test_make_aligner_invalid_gap_scores():
    with pytest.raises(ValueError):
        # open_gap_score > extend_gap_score should error
        make_aligner(open_gap_score=0, extend_gap_score=-1)


@pytest.mark.parametrize(
    "seq1, seq2, blocks1, blocks2, exp1, exp2",
    [
        # perfect match
        ("ACGT", "ACGT", [(0, 4)], [(0, 4)], "ACGT", "ACGT"),
        # case: simple one-gap in seq2
        ("ACGT", "AGGT", [(0, 1), (1, 4)], [(0, 1), (1, 4)], "ACGT", "AGGT"),
        # case: deletion in seq2 → a gap
        ("ACGT", "AGT", [(0, 1), (2, 4)], [(0, 1), (1, 3)], "ACGT", "A-GT"),
    ],
)
def test_apply_alignment_gaps_basic(
    seq1, seq2, blocks1, blocks2, exp1, exp2
):
    out1, out2 = apply_alignment_gaps(seq1, seq2, blocks1, blocks2)
    assert out1 == exp1
    assert out2 == exp2


@pytest.mark.parametrize(
    "seq1, seq2, blocks1, blocks2",
    [
        # mismatched block list lengths
        ("ACGT", "ACGT", [(0, 2)], [(0, 2), (2, 4)]),
        # zero-length block
        ("AC", "AC", [(0, 0)], [(0, 1)]),
        # out-of-bounds
        ("AC", "AC", [(0, 3)], [(0, 1)]),
        # overlapping / unsorted
        ("AC", "AC", [(1, 2), (0, 1)], [(1, 2), (0, 1)]),
    ],
)
def test_apply_alignment_gaps_errors(seq1, seq2, blocks1, blocks2):
    with pytest.raises(ValueError):
        apply_alignment_gaps(seq1, seq2, blocks1, blocks2)


# ─── Tests for SequenceCoordinateMapper ────────────────────────────────────

@pytest.fixture
def mapper():
    return SequenceCoordinateMapper()


def test_find_exon_alignment_blocks_global(mapper):
    # when exon_order is None, does a single global block
    seq = "ATGC"
    b_sp, b_tx = mapper._find_exon_alignment_blocks(seq, seq, None, None)
    assert b_sp == [(0, 4)]
    assert b_tx == [(0, 4)]


def test_find_exon_alignment_blocks_per_exon(mapper):
    # exact‐match branch: exon seq found by simple find()
    spliced = "AAA"
    transcript = "XXXAAAYYY"
    exon_defs = {1: {"exon_number": 1, "start": 1, "end": 3, "sequence": "AAA"}}
    b_sp, b_tx = mapper._find_exon_alignment_blocks(
        spliced, transcript, [1], exon_defs
    )
    assert b_sp == [(0, 3)]
    assert b_tx == [(3, 6)]


def test_extract_base_mappings_simple_run(mapper):
    sp = "ATGC"
    tx = "ATGC"
    # raw coords with positions = indices
    raw_g = [RawBase("r", i, b) for i, b in enumerate(sp)]
    raw_t = [RawBase("r", i, b) for i, b in enumerate(tx)]
    blocks = [(0, 4)]
    # with min_block_length=1, all four should map
    result = mapper._extract_base_mappings(
        sp, tx, raw_g, raw_t, blocks, blocks, min_block_length=1
    )
    assert len(result) == 4
    for idx, m in enumerate(result):
        assert isinstance(m, BaseCoordinateMapping)
        assert m.genomic_base == sp[idx]
        assert m.transcript_base == tx[idx]
        assert m.genomic_position == m.transcript_position == idx


def test_extract_base_mappings_respects_min_block(mapper):
    sp = "ATGC"
    tx = "ATGC"
    raw_g = [RawBase("r", i, b) for i, b in enumerate(sp)]
    raw_t = [RawBase("r", i, b) for i, b in enumerate(tx)]
    blocks = [(0, 4)]
    # min_block_length > run length → no mappings
    result = mapper._extract_base_mappings(
        sp, tx, raw_g, raw_t, blocks, blocks, min_block_length=5
    )
    assert result == []


def test_generate_codon_mappings_triples(mapper):
    # build 6 bases → two codons: ATG (M), CGA (R)
    seq = "ATGCGA"
    base_maps = [
        BaseCoordinateMapping(i, i, i, b, b)
        for i, b in enumerate(seq)
    ]
    cod_map, dna2aa, cdna2aa = mapper._generate_codon_mappings(
        base_maps, seq_region="r", transcript_id="tx"
    )
    # two codons
    assert list(cod_map.keys()) == [1, 2]
    assert cod_map[1].codon == "ATG"
    assert cod_map[1].amino_acid == "M"
    assert cod_map[2].codon == "CGA"
    assert cod_map[2].amino_acid == "R"
    # dna2aa and cdna2aa map each of the three positions → correct AA
    for pos in (0, 1, 2):
        assert dna2aa[pos] == "M"
        assert cdna2aa[pos] == "M"
    for pos in (3, 4, 5):
        assert dna2aa[pos] == "R"
        assert cdna2aa[pos] == "R"


def test_map_transcripts_simple_case(mapper):
    meta = {"seq_region_accession": "chrTest", "start": 1, "strand": 1}
    full_genomic = "ATGCGA"
    # one transcript with one exon covering the entire region
    exon_defs = {
        "tx1": [
            {"exon_number": 1, "start": 1, "end": 6, "sequence": full_genomic}
        ]
    }
    tx_seqs = {"tx1": full_genomic}
    exon_orders = {"tx1": [1]}

    results = mapper.map_transcripts(
        genome_metadata=meta,
        full_genomic_sequence=full_genomic,
        exon_definitions_by_transcript=exon_defs,
        transcript_sequences=tx_seqs,
        exon_orders=exon_orders,
        min_block_length=1,
    )
    # only tx1 should appear
    assert set(results.keys()) == {"tx1"}
    res: TranscriptMappingResult = results["tx1"]

    # base‐level mapping: 6 bases
    assert len(res.transcript_to_genomic) == 6
    # 1-based cdna<->genomic
    assert res.cdna_to_dna_map == {i: i for i in range(1, 7)}
    assert res.dna_to_cdna_map == {i: i for i in range(1, 7)}
    # codon map and protein map
    assert res.codon_map[1].amino_acid == "M"  # ATG
    assert res.codon_map[2].amino_acid == "R"  # CGA
    # gapped sequences should be ungapped
    assert res.gapped_full_genome_sequence[1] == full_genomic
    assert res.gapped_transcript_sequence[1] == full_genomic


def test_map_transcripts_skips_missing_defs(mapper):
    meta = {"seq_region_accession": "chrTest", "start": 1, "strand": 1}
    full_genomic = "ATGCGA"
    # only define exons for tx1
    exon_defs = {
        "tx1": [
            {"exon_number": 1, "start": 1, "end": 6, "sequence": full_genomic}
        ]
    }
    tx_seqs = {"tx1": full_genomic, "tx2": full_genomic}
    exon_orders = {"tx1": [1], "tx2": [1]}

    results = mapper.map_transcripts(
        genome_metadata=meta,
        full_genomic_sequence=full_genomic,
        exon_definitions_by_transcript=exon_defs,
        transcript_sequences=tx_seqs,
        exon_orders=exon_orders,
        min_block_length=1,
    )
    assert set(results.keys()) == {"tx1"}


def test_map_transcripts_exon_length_mismatch_raises(mapper):
    meta = {"seq_region_accession": "chrTest", "start": 1, "strand": 1}
    full_genomic = "ATGCGA"
    # exon end-start+1 != len(sequence) → should raise in build_raw_genome_coords
    exon_defs = {
        "tx1": [
            {"exon_number": 1, "start": 1, "end": 5, "sequence": full_genomic}
        ]
    }
    tx_seqs = {"tx1": full_genomic}
    exon_orders = {"tx1": [1]}

    with pytest.raises(ValueError):
        mapper.map_transcripts(
            genome_metadata=meta,
            full_genomic_sequence=full_genomic,
            exon_definitions_by_transcript=exon_defs,
            transcript_sequences=tx_seqs,
            exon_orders=exon_orders,
            min_block_length=1,
        )
