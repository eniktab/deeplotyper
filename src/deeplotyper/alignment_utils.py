from __future__ import annotations

from functools import cache
from collections.abc import Sequence
from typing import TypeAlias, Tuple

from Bio.Align import PairwiseAligner

__all__ = ("make_aligner", "apply_alignment_gaps")

# Constants & Types
Block: TypeAlias = Tuple[int, int]
GAP_CHAR = "-"

DEFAULT_MATCH_SCORE = 5.0
DEFAULT_MISMATCH_SCORE = -4.0
DEFAULT_OPEN_GAP_SCORE = -10.0
DEFAULT_EXTEND_GAP_SCORE = -1.0


# Exceptions
class AlignmentError(ValueError):
    """Base class for alignment‐related errors."""


class BlockValidationError(AlignmentError):
    """Raised when an aligned‐block interval is invalid."""


@cache
def make_aligner(
    mode: str = "local",
    match_score: float = DEFAULT_MATCH_SCORE,
    mismatch_score: float = DEFAULT_MISMATCH_SCORE,
    open_gap_score: float = DEFAULT_OPEN_GAP_SCORE,
    extend_gap_score: float = DEFAULT_EXTEND_GAP_SCORE,
) -> PairwiseAligner:
    """Return a cached PairwiseAligner configured with custom scoring.

    Args:
        mode (str): Alignment mode ('local', 'global', 'semiglobal').
        match_score (float): Score for a match.
        mismatch_score (float): Score for a mismatch.
        open_gap_score (float): Penalty for opening a gap (must be <= extend_gap_score).
        extend_gap_score (float): Penalty for extending a gap.

    Returns:
        PairwiseAligner: Configured Biopython aligner.

    Raises:
        ValueError: If open_gap_score > extend_gap_score.
    """
    if open_gap_score > extend_gap_score:
        raise ValueError(
            f"open_gap_score ({open_gap_score}) must be <= "
            f"extend_gap_score ({extend_gap_score})"
        )

    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    return aligner


def _validate_blocks(blocks: Sequence[Block], seq_len: int, name: str) -> None:
    """Ensure blocks are sorted, non-overlapping, non-empty, and within bounds.

    Args:
        blocks (Sequence[Block]): Sequence of (start, end) pairs.
        seq_len (int): Length of the sequence being validated.
        name (str): Identifier for error messages.

    Raises:
        BlockValidationError: If any block is zero-length, out of order, overlaps,
            or lies outside 0..seq_len.
    """
    prev_end = 0
    for i, (start, end) in enumerate(blocks):
        if start == end:
            raise BlockValidationError(f"{name}[{i}] has zero length")
        if start < prev_end:
            raise BlockValidationError(f"{name}[{i}] overlaps or is out of order")
        if not (0 <= start <= end <= seq_len):
            raise BlockValidationError(
                f"{name}[{i}] = ({start}, {end}) outside 0..{seq_len}"
            )
        prev_end = end


def apply_alignment_gaps(
    seq1: str, seq2: str, blocks1: Sequence[Block], blocks2: Sequence[Block]
) -> tuple[str, str]:
    """Insert gaps so aligned blocks line up between two sequences.

    Args:
        seq1 (str): First (reference) sequence.
        seq2 (str): Second sequence.
        blocks1 (Sequence[Block]): Aligned‐block intervals for seq1.
        blocks2 (Sequence[Block]): Aligned‐block intervals for seq2.

    Returns:
        tuple[str, str]: Two gapped sequences in which each aligned block
            appears at the same indices in both outputs.

    Raises:
        AlignmentError: If blocks1 and blocks2 differ in length.
        BlockValidationError: If any block is invalid (zero-length,
            overlapping, or out of bounds).
    """
    if len(blocks1) != len(blocks2):
        raise AlignmentError(
            f"block count mismatch: len(blocks1)={len(blocks1)} "
            f"!= len(blocks2)={len(blocks2)}"
        )

    _validate_blocks(blocks1, len(seq1), "blocks1")
    _validate_blocks(blocks2, len(seq2), "blocks2")

    parts1: list[str] = []
    parts2: list[str] = []
    cursor1 = cursor2 = 0

    for (start1, end1), (start2, end2) in zip(blocks1, blocks2):
        # Unaligned leading segment
        if start1 > cursor1:
            parts1.append(seq1[cursor1:start1])
            parts2.append(GAP_CHAR * (start1 - cursor1))
        if start2 > cursor2:
            parts1.append(GAP_CHAR * (start2 - cursor2))
            parts2.append(seq2[cursor2:start2])

        # Aligned block
        parts1.append(seq1[start1:end1])
        parts2.append(seq2[start2:end2])
        cursor1, cursor2 = end1, end2

    # Trailing tails
    tail1_len = len(seq1) - cursor1
    tail2_len = len(seq2) - cursor2

    if tail1_len > 0:
        parts1.append(seq1[cursor1:])
        parts2.append(GAP_CHAR * tail1_len)
    if tail2_len > 0:
        parts1.append(GAP_CHAR * tail2_len)
        parts2.append(seq2[cursor2:])

    return "".join(parts1), "".join(parts2)
