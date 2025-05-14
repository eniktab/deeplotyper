from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Optional

from .data_models import RawBase, ExonDef

def build_linear_coords(
    sequence: str,
    region: str,
    start_pos: int,
    strand: int,
) -> list[RawBase]:
    """Map a linear sequence to a list of RawBase objects.

    Args:
        sequence: The nucleotide sequence.
        region: Name of the region (will be assigned to every RawBase.region).
        start_pos: Genomic or transcript position of the first base.
        strand: +1 for forward strand, -1 for reverse strand.

    Returns:
        A list of RawBase, in ascending coordinate order if strand=+1,
        or descending if strand=-1.

    Raises:
        ValueError: If `strand` is not +1 or -1.
    """
    if strand not in (1, -1):
        raise ValueError(f"strand must be +1 or -1, got {strand}")
    seq = sequence if strand == 1 else sequence[::-1]
    return [
        RawBase(region=region, position=start_pos + i, base=base)
        for i, base in enumerate(seq)
    ]


def _build_exon_coords(
    exon_order: Sequence[int],
    exon_defs: Mapping[int, ExonDef],
    region: str,
    strand: int,
) -> list[RawBase]:
    """Helper to assemble RawBase from ordered exons.

    Args:
        exon_order: List of exon IDs in the desired order.
        exon_defs: Mapping from exon ID to an ExonDef dict-like with keys
            'sequence', 'start', and 'end'.
        region: Region name to assign to each RawBase.
        strand: +1 or -1 orientation.

    Returns:
        A list of RawBase spanning all exons in `exon_order`.

    Raises:
        ValueError: If any exon’s sequence length doesn’t match its coordinates.
    """
    coords: list[RawBase] = []
    for ex_id in exon_order:
        exon = exon_defs[ex_id]
        seq = exon["sequence"]
        start, end = exon["start"], exon["end"]
        positions = list(range(start, end + 1))
        if len(seq) != len(positions):
            raise ValueError(
                f"Exon {ex_id} sequence length {len(seq)} "
                f"!= coordinate span {len(positions)}"
            )
        if strand < 0:
            positions.reverse()
        coords.extend(
            RawBase(region=region, position=pos, base=b)
            for pos, b in zip(positions, seq)
        )
    return coords


def build_raw_genome_coords(
    spliced: str,
    region: str,
    strand: int,
    start_offset: int,
    exon_order: Optional[Sequence[int]] = None,
    exon_defs: Optional[Mapping[int, ExonDef]] = None,
) -> list[RawBase]:
    """Generate RawBase list for a spliced genome region.

    If both `exon_order` and `exon_defs` are provided, builds per–exon
    coordinates; otherwise treats `spliced` as a contiguous run.

    Args:
        spliced: The concatenated (spliced) sequence.
        region: Name of the sequence region.
        strand: +1 for forward or -1 for reverse.
        start_offset: Starting coordinate if falling back to linear mapping.
        exon_order: Optional ordered exon IDs.
        exon_defs: Optional mapping from exon ID to ExonDef.

    Returns:
        A list of RawBase covering either each exon or the full spliced run.

    Raises:
        ValueError: If only one of `exon_order` or `exon_defs` is provided.
    """
    if (exon_order is not None) ^ (exon_defs is not None):
        raise ValueError(
            "Must provide both exon_order and exon_defs, or neither."
        )

    if exon_order is not None:
        return _build_exon_coords(exon_order, exon_defs, region, strand)

    # Fallback: simple linear mapping
    return [
        RawBase(region=region, position=start_offset + i, base=b)
        for i, b in enumerate(spliced)
    ]


def build_raw_transcript_coords(
    transcript_seq: str,
    transcript_name: str,
    transcript_start: int,
) -> list[RawBase]:
    """Generate RawBase list for a transcript (always forward strand).

    Args:
        transcript_seq: The transcript nucleotide sequence.
        transcript_name: Name to assign to each RawBase.region.
        transcript_start: Starting transcript coordinate.

    Returns:
        A list of RawBase mapping the transcript.
    """
    return build_linear_coords(
        sequence=transcript_seq,
        region=transcript_name,
        start_pos=transcript_start,
        strand=1,
    )
