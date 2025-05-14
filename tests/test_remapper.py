import pytest

from deeplotyper import make_aligner, apply_alignment_gaps
from deeplotyper import (
    build_linear_coords,
    build_raw_genome_coords,
    build_raw_transcript_coords,
)
from deeplotyper.orf_utils import find_orfs, get_longest_orf
from deeplotyper.mapper import SequenceCoordinateMapper
from deeplotyper import HaplotypeRemapper
from deeplotyper.data_models import HaplotypeEvent, RawBase


# ─── Sample data ────────────────────────────────────────────────────────────────

genomic_seq = "ATG" + "TT" + "AAA" + "T" + "TAG"  # "ATGTTAAATTAG"
genomic_info = {
    "seq_region_accession": "chr1",
    "start": 80,
    "strand": 1,
}

exon_info_by_tx = {
    "tx1": [
        {"exon_number": 1, "start": 80, "end": 82, "strand": 1, "sequence": "ATG"},
        {"exon_number": 3, "start": 89, "end": 91, "strand": 1, "sequence": "TAG"},
    ],
    "tx2": [
        {"exon_number": 1, "start": 80, "end": 82, "strand": 1, "sequence": "ATG"},
        {"exon_number": 2, "start": 85, "end": 87, "strand": 1, "sequence": "AAA"},
        {"exon_number": 3, "start": 89, "end": 91, "strand": 1, "sequence": "TAG"},
    ],
}

exon_orders = {
    "tx1": [1, 3],
    "tx2": [1, 2, 3],
}

transcript_seqs = {
    "tx1": "ATGTAG",
    "tx2": "ATGAAATAG",
}

@pytest.fixture(scope="module")
def mapping_results():
    return SequenceCoordinateMapper().map_transcripts(
        genome_metadata=genomic_info,
        full_genomic_sequence=genomic_seq,
        exon_definitions_by_transcript=exon_info_by_tx,
        transcript_sequences=transcript_seqs,
        exon_orders=exon_orders,
        min_block_length=3,
    )


# ─── Existing mapping tests ─────────────────────────────────────────────────────

def test_mapping_metadata_and_roundtrip_maps(mapping_results):
    for tx_id, tx in mapping_results.items():
        # 1) metadata
        assert tx.seq_region == "chr1"
        assert tx.offset == 80

        # 2) cdna_to_dna_map and dna_to_cdna_map are exact inverses
        inv = {v: k for k, v in tx.cdna_to_dna_map.items()}
        assert tx.dna_to_cdna_map == inv

        # 3) exon_to_genomic covers exactly the same positions as the cDNA map
        exon_positions = sorted(
            m["genomic_index"]
            for mappings in tx.exon_to_genomic.values()
            for m in mappings
        )
        cdna_positions = sorted(tx.cdna_to_dna_map.values())
        assert exon_positions == cdna_positions

        # 4) codon_map codons line up with the transcript sequence
        seq = transcript_seqs[tx_id]
        expected_codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
        actual_codons = [cm.codon for _, cm in sorted(tx.codon_map.items())]
        assert actual_codons == expected_codons


@pytest.mark.parametrize(
    "events, expected_genome, expected_tx_seqs",
    [
        ((HaplotypeEvent(2, "G", "A"),), "ATATTAAATTAG", {
            "tx1": "ATATAG", "tx2": "ATAAAATAG"
        }),
        ((HaplotypeEvent(3, "", "C"),), "ATGTCTAAATTAG", {
            "tx1": "ATGTAG", "tx2": "ATGAAATAG"
        }),
        ((HaplotypeEvent(5, "A", ""),), "ATGTTAATTAG", {
            "tx1": "ATGTAG", "tx2": "ATGAATAG"
        }),
        (
            (HaplotypeEvent(1, "T", "C"), HaplotypeEvent(7, "A", "")),
            "ACGTTAATTAG", {
                "tx1": "ACGTAG", "tx2": "ACGAATAG"
            }
        ),
        ((HaplotypeEvent(6, "", "GGGGG"),), "ATGTTAGGGGGAATTAG", {
            "tx1": "ATGTAG", "tx2": "ATGAAGGGGGATAG"
        }),
        ((HaplotypeEvent(4, "TAAAT", ""),), "ATGTTAG", {
            "tx1": "ATGTAG", "tx2": "ATGTAG"
        }),
    ]
)
def test_apply_haplotypes_full(mapping_results, events, expected_genome, expected_tx_seqs):
    remapper = HaplotypeRemapper(genomic_seq, mapping_results)
    mutated = remapper.apply_haplotypes({events: ("sample",)})
    out = mutated[frozenset(events)]

    # A) required keys
    assert set(out.keys()) == {
        "mutated_genome",
        "genome_protein",
        "orig_to_mutated_index",
        "mutated_index_to_orig",
        "transcripts"
    }

    # B) mutated genome
    assert out["mutated_genome"] == expected_genome

    # C) coordinate‐map invertibility
    o2m = out["orig_to_mutated_index"]
    m2o = out["mutated_index_to_orig"]
    for orig_pos, mut_idx in o2m.items():
        assert m2o[mut_idx] == orig_pos

    # D) all transcripts present
    assert set(out["transcripts"].keys()) == set(mapping_results.keys())

    # E) per‐transcript cDNA and non‐empty genome_protein
    assert isinstance(out["genome_protein"], str) and out["genome_protein"]
    for tx_id, exp_cdna in expected_tx_seqs.items():
        assert out["transcripts"][tx_id].cdna_sequence == exp_cdna


# ─── New edge‐case tests ─────────────────────────────────────────────────────────

def test_make_aligner_caching_and_modes():
    a1 = make_aligner()
    a2 = make_aligner()
    assert a1 is a2
    a3 = make_aligner("global")
    a4 = make_aligner("global")
    assert a3 is a4
    assert a1.mode == "local"
    assert a3.mode == "global"


def test_apply_alignment_gaps_simple_internal_gap():
    seq1 = "ATCG"
    seq2 = "ACG"
    blocks1 = [(0, 1), (2, 4)]
    blocks2 = [(0, 1), (1, 3)]
    g1, g2 = apply_alignment_gaps(seq1, seq2, blocks1, blocks2)
    assert g1 == "ATCG"
    assert g2 == "A-CG"


def test_build_linear_and_transcript_coords():
    seq = "ABC"
    coords_p = build_linear_coords(seq, "reg", 10, strand=1)
    expected_p = [
        RawBase("reg", 10, "A"),
        RawBase("reg", 11, "B"),
        RawBase("reg", 12, "C"),
    ]
    assert coords_p == expected_p

    coords_n = build_linear_coords(seq, "reg", 10, strand=-1)
    expected_n = [
        RawBase("reg", 10, "C"),
        RawBase("reg", 11, "B"),
        RawBase("reg", 12, "A"),
    ]
    assert coords_n == expected_n

    tx_coords = build_raw_transcript_coords("XYZ", "tx", 5)
    expected_tx = [
        RawBase("tx", 5, "X"),
        RawBase("tx", 6, "Y"),
        RawBase("tx", 7, "Z"),
    ]
    assert tx_coords == expected_tx


def test_build_raw_genome_coords_fallback_and_exon_mismatch():
    fallback = build_raw_genome_coords("ATG", "chrX", 1, 100, None, None)
    assert fallback == [
        RawBase("chrX", 100, "A"),
        RawBase("chrX", 101, "T"),
        RawBase("chrX", 102, "G"),
    ]

    bad_defs = {1: {"exon_number": 1, "start": 1, "end": 4, "sequence": "ATG"}}
    with pytest.raises(ValueError):
        build_raw_genome_coords("ATG", "chrX", 1, 1, [1], bad_defs)


def test_find_orfs_and_get_longest_orf():
    seq = "CCCATGAAATGATAGGGTGA"
    orfs = find_orfs(seq)
    assert any(start == 3 and seq[start:end] == seq[3:end] for start, end, _ in orfs)
    assert any(start == 8 and seq[start:end] == seq[8:end] for start, end, _ in orfs)

    longest = get_longest_orf(seq)
    assert longest == max(orfs, key=lambda x: len(x[2]))

    seq2 = "CCCCCCC"
    assert find_orfs(seq2) == []
    assert get_longest_orf(seq2) == (0, len(seq2), seq2)


def test_skips_undefined_transcripts():
    mapper = SequenceCoordinateMapper()
    results = mapper.map_transcripts(
        genome_metadata=genomic_info,
        full_genomic_sequence=genomic_seq,
        exon_definitions_by_transcript={"tx1": exon_info_by_tx["tx1"]},
        transcript_sequences={"tx1": "ATGTAG", "tx3": "AAAA"},
        exon_orders={"tx1": [1, 3]},
        min_block_length=3,
    )
    assert "tx1" in results
    assert "tx3" not in results


def test_mapping_negative_strand_behavior():
    mapper = SequenceCoordinateMapper()
    gm = "ATG"
    gm_info = {"seq_region_accession": "chr1", "start": 1, "strand": -1}
    defs = {"tx": [{"exon_number": 1, "start": 1, "end": 3, "sequence": "ATG"}]}
    orders = {"tx": [1]}
    seqs = {"tx": "ATG"}

    res = mapper.map_transcripts(
        genome_metadata=gm_info,
        full_genomic_sequence=gm,
        exon_definitions_by_transcript=defs,
        transcript_sequences=seqs,
        exon_orders=orders,
        min_block_length=1,
    )
    tx = res["tx"]
    pos_list = [m["genomic_index"] for m in tx.transcript_to_genomic]
    assert pos_list == [3, 2, 1]


def test_apply_haplotypes_empty_and_bad_ref(mapping_results):
    remapper = HaplotypeRemapper(genomic_seq, mapping_results)
    assert remapper.apply_haplotypes({}) == {}

    bad = {(HaplotypeEvent(0, "C", "G"),): ("sample",)}
    with pytest.raises(AssertionError):
        remapper.apply_haplotypes(bad)


# ─── Additional edge‐case tests for apply_haplotypes ───────────────────────────

def test_apply_haplotypes_insertion_at_start(mapping_results):
    remapper = HaplotypeRemapper(genomic_seq, mapping_results)
    alt = "GG"
    ev = (HaplotypeEvent(0, "", alt),)
    out = remapper.apply_haplotypes({ev: ("ins_start",)})[frozenset(ev)]

    assert out["mutated_genome"] == alt + genomic_seq
    for tx_id, tx in out["transcripts"].items():
        orig = transcript_seqs[tx_id]
        expected = orig[0] + alt + orig[1:]
        assert tx.cdna_sequence == expected


def test_apply_haplotypes_deletion_at_end(mapping_results):
    remapper = HaplotypeRemapper(genomic_seq, mapping_results)
    pos = len(genomic_seq) - 1
    ref = genomic_seq[pos]
    ev = (HaplotypeEvent(pos, ref, ""),)
    out = remapper.apply_haplotypes({ev: ("del_end",)})[frozenset(ev)]

    assert out["mutated_genome"] == genomic_seq[:-1]
    for tx_id, tx in out["transcripts"].items():
        assert tx.cdna_sequence == transcript_seqs[tx_id][:-1]


def test_apply_haplotypes_multibase_substitution(mapping_results):
    remapper = HaplotypeRemapper(genomic_seq, mapping_results)
    pos = 2
    ref = genomic_seq[pos : pos + 3]
    alt = "CCC"
    ev = (HaplotypeEvent(pos, ref, alt),)
    out = remapper.apply_haplotypes({ev: ("multi_sub",)})[frozenset(ev)]

    expected_genome = genomic_seq[:pos] + alt + genomic_seq[pos + 3 :]
    assert out["mutated_genome"] == expected_genome

    # transcripts pick up exactly those three bases in exon1
    assert out["transcripts"]["tx1"].cdna_sequence == "ATCTAG"
    assert out["transcripts"]["tx2"].cdna_sequence == "ATCAAATAG"


def test_apply_haplotypes_multiple_substitutions(mapping_results):
    remapper = HaplotypeRemapper(genomic_seq, mapping_results)
    evs = (
        HaplotypeEvent(1, genomic_seq[1], "C"),
        HaplotypeEvent(4, genomic_seq[4], "G"),
    )
    out = remapper.apply_haplotypes({evs: ("dual_snp",)})[frozenset(evs)]
    mg = out["mutated_genome"]

    # genome has both SNPs
    assert mg == "ACGTGAAATTAG"

    # transcripts accumulate only the overlapping SNP
    for tx_id, tx in out["transcripts"].items():
        # build expected by applying only pos=1 change
        seq = list(transcript_seqs[tx_id])
        dna2cdna = mapping_results[tx_id].dna_to_cdna_map
        gp1 = mapping_results[tx_id].offset + 1
        if gp1 in dna2cdna:
            idx = dna2cdna[gp1] - 1
            seq[idx] = "C"
        assert "".join(seq) == tx.cdna_sequence


def test_apply_haplotypes_deletion_in_intron(mapping_results):
    remapper = HaplotypeRemapper(genomic_seq, mapping_results)
    ev = (HaplotypeEvent(4, genomic_seq[4], ""),)
    out = remapper.apply_haplotypes({ev: ("del_intron",)})[frozenset(ev)]

    assert out["mutated_genome"] == genomic_seq[:4] + genomic_seq[5:]
    for tx_id, tx in out["transcripts"].items():
        assert tx.cdna_sequence == transcript_seqs[tx_id]


def test_apply_haplotypes_insertion_in_intron(mapping_results):
    remapper = HaplotypeRemapper(genomic_seq, mapping_results)
    alt = "GGG"
    ev = (HaplotypeEvent(4, "", alt),)
    out = remapper.apply_haplotypes({ev: ("ins_intron",)})[frozenset(ev)]

    assert out["mutated_genome"] == genomic_seq[:4] + alt + genomic_seq[4:]
    for tx_id, tx in out["transcripts"].items():
        assert tx.cdna_sequence == transcript_seqs[tx_id]


def test_apply_haplotypes_full_deletion(mapping_results):
    remapper = HaplotypeRemapper(genomic_seq, mapping_results)
    ev = (HaplotypeEvent(0, genomic_seq, ""),)
    out = remapper.apply_haplotypes({ev: ("del_all",)})[frozenset(ev)]

    assert out["mutated_genome"] == ""
    for tx_id, tx in out["transcripts"].items():
        assert tx.cdna_sequence == ""


def test_apply_haplotypes_overlapping_events_raise(mapping_results):
    remapper = HaplotypeRemapper(genomic_seq, mapping_results)
    ev = (
        HaplotypeEvent(2, genomic_seq[2:5], "A"),
        HaplotypeEvent(3, genomic_seq[3], "C"),
    )
    # overlapping multi‐base deletion is applied sequentially
    out = remapper.apply_haplotypes({ev: ("overlap",)})[frozenset(ev)]

    # first drop "GTT" → "ATAAAATTAG", then at new index 1 replace "T" with "C"
    assert out["mutated_genome"] == "ACAAAATTAG"


def test_apply_haplotypes_invalid_positions_and_mismatch(mapping_results):
    remapper = HaplotypeRemapper(genomic_seq, mapping_results)

    # negative-position insertion prepends
    ev_neg = (HaplotypeEvent(-1, "", "A"),)
    out_neg = remapper.apply_haplotypes({ev_neg: ("neg_ins",)})[frozenset(ev_neg)]
    assert out_neg["mutated_genome"] == "A" + genomic_seq

    # substitution past end => error
    with pytest.raises(AssertionError):
        remapper.apply_haplotypes({(HaplotypeEvent(len(genomic_seq), "A", "T"),): ("oob_snp",)})

    # insertion past end => append
    ev_oob = (HaplotypeEvent(len(genomic_seq) + 1, "", "A"),)
    out_oob = remapper.apply_haplotypes({ev_oob: ("oob_ins",)})[frozenset(ev_oob)]
    assert out_oob["mutated_genome"] == genomic_seq + "A"

    # ref‐allele mismatch => error
    with pytest.raises(AssertionError):
        remapper.apply_haplotypes({(HaplotypeEvent(2, "AAA", "CCC"),): ("bad_ref",)})
