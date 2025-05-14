from __future__ import annotations

from dataclasses import dataclass, field
from typing import TypedDict, Mapping, List, Tuple, TypeAlias

# Type aliases for clarity
Position: TypeAlias = int
Base: TypeAlias = str
CodonIndex: TypeAlias = int
Region: TypeAlias = str

class ExonDef(TypedDict):
    """
    Definition of an exon segment.
    """
    exon_number: int
    start: Position
    end: Position
    sequence: str

class MappingEntry(TypedDict):
    """
    Mapping between transcript and genomic indices.
    """
    transcript_index: Position
    transcript_base: Base
    genomic_index: Position
    genomic_base: Base

@dataclass(order=True, frozen=True, slots=True)
class BaseCoordinateMapping:
    """
    Maps a single base between global, genomic, and transcript coordinates.

    Attributes:
        global_position (Position): 0-based index in the global coordinate space.
        genomic_position (Position): 0-based index in the genomic sequence.
        transcript_position (Position): 0-based index in the transcript.
        genomic_base (Base): Nucleotide at genomic_position.
        transcript_base (Base): Nucleotide at transcript_position.
    """
    global_position: Position
    genomic_position: Position
    transcript_position: Position
    genomic_base: Base
    transcript_base: Base

@dataclass(frozen=True, slots=True)
class CodonCoordinateMapping:
    """
    Maps a codon and its amino acid across genomic and transcript coordinates.

    Attributes:
        codon_index (CodonIndex): Index of the codon in the transcript (0-based).
        codon (str): Three-character codon sequence.
        amino_acid (str): Single-letter amino acid code.
        genomic_positions (List[Tuple[Region, Position]]):
            List of (region, position) for each base in the codon on the genome.
        transcript_positions (List[Tuple[Region, Position]]):
            List of (region, position) for each base in the codon on the transcript.
    """
    codon_index: CodonIndex
    codon: str
    amino_acid: str
    genomic_positions: List[Tuple[Region, Position]]
    transcript_positions: List[Tuple[Region, Position]]

@dataclass
class SequenceMappingResult:
    """
    Result of mapping an entire transcript against the genome at base resolution.

    Attributes:
        base_mappings (List[BaseCoordinateMapping]):
            Ordered list of base-wise mappings.
        base_mappings_by_genome (List[BaseCoordinateMapping]):
            Same mappings sorted by genomic_position.
        gapped_full_genome_sequence (Mapping[Position, Base]):
            Genome sequence with alignment gaps inserted.
        gapped_transcript_sequence (Mapping[Position, Base]):
            Transcript sequence with alignment gaps inserted.
        gapped_alignment_map (Mapping[Position, List[Tuple[Position, Position]]]):
            For each gapped genome index, list of (gapped_transcript_index, length) pairs.
    """
    base_mappings: List[BaseCoordinateMapping] = field(default_factory=list)
    base_mappings_by_genome: List[BaseCoordinateMapping] = field(default_factory=list)
    gapped_full_genome_sequence: Mapping[Position, Base] = field(default_factory=dict)
    gapped_transcript_sequence: Mapping[Position, Base] = field(default_factory=dict)
    gapped_alignment_map: Mapping[Position, List[Tuple[Position, Position]]] = field(default_factory=dict)

@dataclass
class TranscriptMappingResult:
    """
    Comprehensive mapping results for a single transcript.

    Attributes:
        transcript_to_genomic (List[MappingEntry]):
            Base-level mappings from transcript indices → genomic indices.
        cdna_to_dna_map (Mapping[Position, Position]):
            cDNA index → genomic DNA index.
        dna_to_cdna_map (Mapping[Position, Position]):
            Genomic DNA index → cDNA index.
        exon_to_genomic (Mapping[int, List[MappingEntry]]):
            For each exon number, list of MappingEntry in genomic coords.
        exon_to_transcript (Mapping[int, List[MappingEntry]]):
            For each exon number, list of MappingEntry in transcript coords.
        dna_to_exon_map (Mapping[Position, Tuple[int, Position]]):
            Genomic DNA index → (exon_number, position_within_exon).
        cdna_to_exon_map (Mapping[Position, Tuple[int, Position]]):
            cDNA index → (exon_number, position_within_exon).
        codon_map (Mapping[CodonIndex, CodonCoordinateMapping]):
            Mapping of each codon index to its coordinate mapping.
        dna_to_protein_map (Mapping[Position, str]):
            Genomic DNA index → single-letter amino acid.
        cdna_to_protein_map (Mapping[Position, str]):
            cDNA index → single-letter amino acid.
        gapped_full_genome_sequence (Mapping[Position, Base]):
            Gapped genome sequence.
        gapped_transcript_sequence (Mapping[Position, Base]):
            Gapped transcript sequence.
        gapped_alignment_map (Mapping[Position, List[Tuple[Position, Position]]]):
            Alignment mapping in gapped coordinates.
        seq_region (str): Name of the genomic region (e.g., chromosome).
        offset (Position): Genomic offset applied to all coordinates.
    """
    transcript_to_genomic: List[MappingEntry] = field(default_factory=list)
    cdna_to_dna_map: Mapping[Position, Position] = field(default_factory=dict)
    dna_to_cdna_map: Mapping[Position, Position] = field(default_factory=dict)
    exon_to_genomic: Mapping[int, List[MappingEntry]] = field(default_factory=dict)
    exon_to_transcript: Mapping[int, List[MappingEntry]] = field(default_factory=dict)
    dna_to_exon_map: Mapping[Position, Tuple[int, Position]] = field(default_factory=dict)
    cdna_to_exon_map: Mapping[Position, Tuple[int, Position]] = field(default_factory=dict)
    codon_map: Mapping[CodonIndex, CodonCoordinateMapping] = field(default_factory=dict)
    dna_to_protein_map: Mapping[Position, str] = field(default_factory=dict)
    cdna_to_protein_map: Mapping[Position, str] = field(default_factory=dict)
    gapped_full_genome_sequence: Mapping[Position, Base] = field(default_factory=dict)
    gapped_transcript_sequence: Mapping[Position, Base] = field(default_factory=dict)
    gapped_alignment_map: Mapping[Position, List[Tuple[Position, Position]]] = field(default_factory=dict)
    seq_region: str = ""
    offset: Position = 0

@dataclass(frozen=True, slots=True)
class RawBase:
    """
    A raw base in a given region before any alignment or mapping.

    Attributes:
        region (Region): Name of the sequence region (e.g., chromosome).
        position (Position): 0-based coordinate in that region.
        base (Base): The nucleotide at that position.
    """
    region: Region
    position: Position
    base: Base

@dataclass(frozen=True, slots=True)
class HaplotypeEvent:
    """
    A single variant event relative to a reference haplotype.

    Attributes:
        pos0 (int): Zero-based offset from the start of the reference window.
        ref_allele (str): Reference allele string at this position.
        alt_seq (str): Alternate allele sequence.
    """
    pos0: int
    ref_allele: str
    alt_seq: str

@dataclass
class NewTranscriptSequences:
    """
    Constructed sequences and maps after applying haplotype events.

    Attributes:
        exon_sequences (Mapping[int, str]):
            Mapping exon_number → exon sequence.
        cdna_sequence (str):
            Full cDNA sequence after events.
        codon_map (Mapping[CodonIndex, CodonCoordinateMapping]):
            Codon coordinate mappings post-mutation.
        dna_to_protein_map (Mapping[Position, str]):
            Genomic DNA index → amino acid.
        cdna_to_protein_map (Mapping[Position, str]):
            cDNA index → amino acid.
    """
    exon_sequences: Mapping[int, str] = field(default_factory=dict)
    cdna_sequence: str = ""
    codon_map: Mapping[CodonIndex, CodonCoordinateMapping] = field(default_factory=dict)
    dna_to_protein_map: Mapping[Position, str] = field(default_factory=dict)
    cdna_to_protein_map: Mapping[Position, str] = field(default_factory=dict)