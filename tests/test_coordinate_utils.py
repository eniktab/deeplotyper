import pytest

from deeplotyper import (
    build_linear_coords,
    build_raw_genome_coords,
    build_raw_transcript_coords,
)
from deeplotyper.data_models import RawBase, ExonDef


def test_build_linear_coords_forward():
    seq = "ATCG"
    region = "chr1"
    start = 10
    result = build_linear_coords(seq, region, start, strand=1)
    expected = [
        RawBase(region="chr1", position=10, base="A"),
        RawBase(region="chr1", position=11, base="T"),
        RawBase(region="chr1", position=12, base="C"),
        RawBase(region="chr1", position=13, base="G"),
    ]
    assert result == expected


def test_build_linear_coords_reverse():
    seq = "ATCG"
    region = "chr2"
    start = 5
    result = build_linear_coords(seq, region, start, strand=-1)
    # Sequence is reversed but positions still ascend
    expected = [
        RawBase(region="chr2", position=5, base="G"),
        RawBase(region="chr2", position=6, base="C"),
        RawBase(region="chr2", position=7, base="T"),
        RawBase(region="chr2", position=8, base="A"),
    ]
    assert result == expected


def test_build_linear_coords_invalid_strand():
    with pytest.raises(ValueError) as exc:
        build_linear_coords("A", "chrX", 0, strand=0)
    assert "strand must be +1 or -1" in str(exc.value)


def test_build_raw_transcript_coords_equivalent():
    # Delegates to build_linear_coords with strand=1
    tx = "GATTACA"
    name = "tx1"
    start = 100
    result = build_raw_transcript_coords(tx, name, start)
    expected = build_linear_coords(tx, name, start, strand=1)
    assert result == expected


def test_build_raw_genome_coords_fallback():
    spliced = "GGCC"
    region = "chrF"
    strand = 1
    offset = 50
    # Neither exon_order nor exon_defs provided → linear fallback
    result = build_raw_genome_coords(
        spliced, region, strand, offset, exon_order=None, exon_defs=None
    )
    expected = [
        RawBase(region="chrF", position=50, base="G"),
        RawBase(region="chrF", position=51, base="G"),
        RawBase(region="chrF", position=52, base="C"),
        RawBase(region="chrF", position=53, base="C"),
    ]
    assert result == expected

    # Empty spliced also yields empty list
    assert build_raw_genome_coords("", region, strand, offset, None, None) == []


@pytest.mark.parametrize("eo,ed", [
    ([1], None),
    (None, {1: {"exon_number": 1, "start": 1, "end": 1, "sequence": "A"}}),
])
def test_build_raw_genome_coords_xor_error(eo, ed):
    with pytest.raises(ValueError) as exc:
        build_raw_genome_coords("A", "chrZ", 1, 0, exon_order=eo, exon_defs=ed)
    assert "Must provide both exon_order and exon_defs" in str(exc.value)


def test_build_raw_genome_coords_single_exon_forward():
    exon_order = [1]
    exon_defs: dict[int, ExonDef] = {
        1: ExonDef(exon_number=1, start=100, end=102, sequence="ACT")
    }
    # spliced is ignored when exons provided
    result = build_raw_genome_coords(
        spliced="", region="chrE", strand=1,
        start_offset=0, exon_order=exon_order, exon_defs=exon_defs
    )
    expected = [
        RawBase(region="chrE", position=100, base="A"),
        RawBase(region="chrE", position=101, base="C"),
        RawBase(region="chrE", position=102, base="T"),
    ]
    assert result == expected


def test_build_raw_genome_coords_single_exon_reverse():
    exon_order = [2]
    exon_defs: dict[int, ExonDef] = {
        2: ExonDef(exon_number=2, start=200, end=202, sequence="GGA")
    }
    result = build_raw_genome_coords(
        spliced="",
        region="chrR",
        strand=-1,
        start_offset=0,
        exon_order=exon_order,
        exon_defs=exon_defs,
    )
    # Positions reversed but sequence order unchanged
    expected = [
        RawBase(region="chrR", position=202, base="G"),
        RawBase(region="chrR", position=201, base="G"),
        RawBase(region="chrR", position=200, base="A"),
    ]
    assert result == expected


def test_build_raw_genome_coords_exon_length_mismatch():
    exon_order = [3]
    exon_defs: dict[int, ExonDef] = {
        3: ExonDef(exon_number=3, start=10, end=12, sequence="AG")
    }
    with pytest.raises(ValueError) as exc:
        build_raw_genome_coords(
            spliced="",
            region="chrM",
            strand=1,
            start_offset=0,
            exon_order=exon_order,
            exon_defs=exon_defs,
        )
    msg = str(exc.value)
    assert "sequence length" in msg and "!= coordinate span" in msg


def test_build_raw_genome_coords_empty_exon_list():
    # Providing empty lists/maps should produce an empty coordinate list
    result = build_raw_genome_coords(
        spliced="XYZ",
        region="chrEmpty",
        strand=1,
        start_offset=0,
        exon_order=[],
        exon_defs={},
    )
    assert result == []
