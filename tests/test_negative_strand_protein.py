"""Decisive minus-strand correctness tests (reverse vs reverse-COMPLEMENT).

These encode a defect the existing suite misses: `test_mapping_negative_strand_behavior`
only checks coordinate ORDER on "ATG" and never verifies a protein, so a strand=-1 path
that REVERSES (rather than leaves intact / reverse-complements) the transcript goes
undetected.

Convention under test (matches the Ensembl fetcher, which returns the SENSE sequence for a
strand=-1 region): the caller supplies SENSE sequences; `strand` governs only the reported
genomic-coordinate DIRECTION. Therefore strand=-1 must yield the SAME cDNA/protein as
strand=+1 for identical sense inputs — differing only in reported genomic positions.

Status: FAILS on the current code (strand=-1 reverses the cDNA). Must PASS after the fix.
"""
from deeplotyper import SequenceCoordinateMapper, HaplotypeRemapper
from deeplotyper.data_models import HaplotypeEvent


def _run(strand, sense, events):
    defs = {"tx": [{"exon_number": 1, "start": 1, "end": len(sense),
                    "strand": strand, "sequence": sense}]}
    mr = SequenceCoordinateMapper().map_transcripts(
        genome_metadata={"seq_region_accession": "g", "start": 1, "strand": strand},
        full_genomic_sequence=sense,
        exon_definitions_by_transcript=defs,
        transcript_sequences={"tx": sense},
        exon_orders={"tx": [1]},
        min_block_length=1,
    )
    rem = HaplotypeRemapper(reference_genome=sense, transcript_results=mr)
    out = rem.apply_haplotypes({events: ("x",)})[frozenset(events)]
    t = out["transcripts"]["tx"]
    protein = "".join(t.codon_map[i].amino_acid for i in sorted(t.codon_map))
    return t.cdna_sequence, protein


def test_negative_strand_matches_positive_strand_reference():
    """Reference (no coding-changing event): strand=-1 must give the same cDNA/protein
    as strand=+1 for identical sense input."""
    sense = "ATGAAATAA"            # M-K-stop  (protein 'MK')
    ev = (HaplotypeEvent(len(sense) - 1, "A", "A"),)   # no-op, to enumerate the haplotype
    cdna_p, prot_p = _run(1, sense, ev)
    cdna_n, prot_n = _run(-1, sense, ev)
    assert cdna_p == sense and prot_p == "MK"
    assert cdna_n == cdna_p, f"strand=-1 cDNA {cdna_n!r} != strand=+1 cDNA {cdna_p!r}"
    assert prot_n == prot_p == "MK"


def test_negative_strand_missense_consequence():
    """A genomic SNV that changes codon 2 (K->E) must yield protein 'ME' on BOTH strands."""
    sense = "ATGAAATAA"
    ev = (HaplotypeEvent(3, "A", "G"),)   # sense pos3(0-based) A>G -> ATGGAATAA -> M-E-stop
    cdna_p, prot_p = _run(1, sense, ev)
    cdna_n, prot_n = _run(-1, sense, ev)
    assert cdna_p == "ATGGAATAA" and prot_p == "ME"
    assert cdna_n == cdna_p, f"strand=-1 cDNA {cdna_n!r} != strand=+1 cDNA {cdna_p!r}"
    assert prot_n == "ME"


def test_negative_strand_insertion_consequence():
    """An in-frame genomic insertion must yield the same cDNA/protein on BOTH strands
    (locks in minus-strand indel reconstruction)."""
    sense = "ATGAAATAA"
    ev = (HaplotypeEvent(2, "", "GGG"),)   # insert a codon -> ATG GGG AAA TAA -> MGK
    cdna_p, prot_p = _run(1, sense, ev)
    cdna_n, prot_n = _run(-1, sense, ev)
    assert cdna_p == "ATGGGGAAATAA" and prot_p == "MGK"
    assert cdna_n == cdna_p, f"strand=-1 cDNA {cdna_n!r} != strand=+1 cDNA {cdna_p!r}"
    assert prot_n == "MGK"
