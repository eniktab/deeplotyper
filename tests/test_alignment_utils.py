import pytest

from deeplotyper.alignment_utils import (apply_alignment_gaps,
                                             AlignmentError,
                                             BlockValidationError)

# ─────────────────────────────────────────────────────────────────────────────
# Valid “flanked exon” cases
# ─────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize(
    "seq1, seq2, blocks1, blocks2, exp1, exp2",
    [
        # 1. Perfect‐match exon in middle of long genomic flanks
        (
            "GGGGAAAATTTTCCCCGGGG",
            "AAAA",
            [(4, 8)],
            [(0, 4)],
            # g1: identical to seq1; g2: only exon, flanked by 4 leading + 12 trailing gaps
            "GGGGAAAATTTTCCCCGGGG",
            "----AAAA------------",
        ),
        # 2. Perfect‐match exon with both leading & trailing flanks
        (
            "NNNNAAAATTTTCCCCMMMM",
            "AAAATTTTCCCC",
            [(4, 16)],  # ← now covers AAAATTTTCCCC in seq1
            [(0, 12)],  # ← covers entire transcript
            "NNNNAAAATTTTCCCCMMMM",
            "----AAAATTTTCCCC----",
        ),
        # 3. Single‐mismatch inside exon, flanked by intronic Ns
        (
            "IIIIAAAACCCCIIII",
            "AAAAGCCC",
            [(4, 12)],
            [(0, 8)],
            "IIIIAAAACCCCIIII",
            "----AAAAGCCC----",
        ),
        # 4. One‐base insertion in transcript, flanked on both sides
        (
            "NNNNATCGTTAANNNN",
            "ATCGGTTAA",
            [(4, 8), (8, 12)],
            [(0, 4), (5, 9)],
            "NNNNATCG-TTAANNNN",
            "----ATCGGTTAA----",
        ),
        # 5. One‐base deletion in transcript, flanked on both sides
        (
            "NNNNATCGTTAANNNN",
            "ACGTTAA",
            [(4, 5), (6, 12)],
            [(0, 1), (1, 7)],
            "NNNNATCGTTAANNNN",
            "----A-CGTTAA----",
        ),
    ],
)
def test_apply_alignment_gaps_flanked_exons(
    seq1, seq2, blocks1, blocks2, exp1, exp2
):
    g1, g2 = apply_alignment_gaps(seq1, seq2, blocks1, blocks2)
    assert g1 == exp1
    assert g2 == exp2


# ─────────────────────────────────────────────────────────────────────────────
# Invalid‐block errors
# ─────────────────────────────────────────────────────────────────────────────
def test_block_count_mismatch_raises():
    with pytest.raises(AlignmentError):
        apply_alignment_gaps("ATCG", "ATCG", [(0, 4)], [(0, 2), (2, 4)])


def test_zero_length_block_raises():
    with pytest.raises(BlockValidationError):
        apply_alignment_gaps("ATCG", "ATCG", [(2, 2)], [(2, 2)])


def test_overlapping_blocks_raises():
    # blocks overlap: (5,10) then (8,12)
    seq = "A" * 20
    with pytest.raises(BlockValidationError):
        apply_alignment_gaps(seq, seq, [(5, 10), (8, 12)], [(5, 10), (8, 12)])


def test_out_of_bounds_block_raises():
    # block end > sequence length
    seq = "A" * 10
    with pytest.raises(BlockValidationError):
        apply_alignment_gaps(seq, seq, [(0, 11)], [(0, 11)])
