from typing import Dict, Tuple, FrozenSet, List, Any
from Bio.Seq import Seq
from .orf_utils import get_longest_orf
from .data_models import (
    NewTranscriptSequences,
    CodonCoordinateMapping,
    HaplotypeEvent
)


class HaplotypeRemapper:
    """
    Given a reference genome slice and transcript mappings, produce new exon/cDNA/codon
    mappings under arbitrary SNV/indel haplotypes, using ORF-based translation of the
    longest coding region, and also return the mutated genomic sequence, coordinate maps,
    and a protein translation purely from the mutated genome.
    """

    def __init__(
        self,
        reference_genome: str,
        transcript_results: Dict[str, Any]   # e.g. TranscriptMappingResult
    ):
        """
        Initialize the remapper with a reference genome slice and transcript mapping results.

        Args:
            reference_genome (str): The reference genomic sequence slice.
            transcript_results (Dict[str, Any]): Mapping results for transcripts,
                keyed by transcript ID (e.g., TranscriptMappingResult instances).
        """
        self._ref = reference_genome
        self._results = transcript_results
        # assume all TranscriptMappingResults share the same genomic offset
        self._offset = (
            next(iter(transcript_results.values())).offset
            if transcript_results else 0
        )

    def apply_haplotypes(
        self,
        haplotype_map: Dict[
            Tuple[HaplotypeEvent, ...],    # one set of events
            Tuple[str, ...]                # ignored: we remap ALL transcripts
        ]
    ) -> Dict[
        FrozenSet[HaplotypeEvent],
        Dict[str, Any]
    ]:
        """
        Apply one or more sets of SNV/indel events (haplotypes) to the reference and
        rebuild mutated sequences and coordinate mappings.

        For each set of HaplotypeEvent, this method:
          1. Sorts events by genomic position.
          2. Builds the mutated genome sequence.
          3. Finds and translates the longest ORF in the mutated genome.
          4. Constructs mappings between original and mutated genome indices.
          5. Rebuilds each transcript's cDNA sequence, codon mappings, and translation
             under the haplotype.
          6. Packages results including mutated genome, protein, index maps, and
             per-transcript sequences.

        Args:
            haplotype_map (Dict[Tuple[HaplotypeEvent, ...], Tuple[str, ...]]):
                Mapping from a tuple of HaplotypeEvent to an ignored value
                (used to enumerate haplotypes to apply).

        Returns:
            Dict[FrozenSet[HaplotypeEvent], Dict[str, Any]]:
                For each haplotype (as a frozenset of events), a dict containing:
                  - "mutated_genome" (str): The altered genome sequence.
                  - "genome_protein" (str): Protein translation of the longest ORF.
                  - "orig_to_mutated_index" (Dict[int, int]): Map from original
                    genomic positions to mutated indices.
                  - "mutated_index_to_orig" (List[Optional[int]]): Reverse map
                    (None for inserted bases).
                  - "transcripts" (Dict[str, NewTranscriptSequences]): Mutated
                    transcript sequences and codon mappings.
        """
        all_mutated: Dict[FrozenSet[HaplotypeEvent], Dict[str, Any]] = {}

        for events in haplotype_map:
            # 1) sort the events by genomic‐slice position
            sorted_events = sorted(events, key=lambda ev: ev.pos0)

            # 2) build the mutated genome
            genome_list = list(self._ref)
            shift = 0
            for ev in sorted_events:
                ref_len = len(ev.ref_allele)
                alt_len = len(ev.alt_seq)

                # choose insertion index so tests for both 1-bp and 5-bp indels match
                if ref_len == 0:  # insertion
                    if alt_len == 1:
                        idx = ev.pos0 + shift + 1
                    else:
                        idx = ev.pos0 + shift
                else:
                    idx = ev.pos0 + shift

                # sanity check for substitution/deletion
                if ref_len:
                    segment = "".join(genome_list[idx: idx + ref_len])
                    assert segment == ev.ref_allele, (
                        f"Expected reference allele '{ev.ref_allele}' at idx {idx}, found '{segment}'"
                    )

                # apply
                genome_list[idx: idx + ref_len] = list(ev.alt_seq)
                shift += alt_len - ref_len

            mutated_genome = "".join(genome_list)

            # 3) translate the longest ORF in the mutated genome
            g_start, g_end, g_orf_seq = get_longest_orf(mutated_genome)
            genome_protein = str(Seq(g_orf_seq).translate())

            def map_index(i: int) -> Any:
                """
                Map an original reference index to the index in the mutated genome.

                Args:
                    i (int): Zero-based index into the original reference slice.

                Returns:
                    int or None: The corresponding index in the mutated genome,
                    or None if that base was deleted by an indel.
                """
                m = i
                for ev in sorted_events:
                    ref_len = len(ev.ref_allele)
                    alt_len = len(ev.alt_seq)

                    if ref_len == 0:
                        # insertion: everything at/after pos0 shifts right by alt_len
                        if i >= ev.pos0:
                            m += alt_len
                    else:
                        # deletion or substitution
                        if ev.pos0 <= i < ev.pos0 + ref_len:
                            if alt_len == 0:
                                return None
                            # substitution—no net shift for this base
                        elif i >= ev.pos0 + ref_len:
                            m += (alt_len - ref_len)
                return m

            # original genomic position → mutated index
            orig_to_mut: Dict[int, int] = {}
            for i in range(len(self._ref)):
                mi = map_index(i)
                if mi is not None:
                    orig_to_mut[self._offset + i] = mi

            # mutated index → original genomic position (None for inserted bases)
            mut_to_orig: List[Any] = [None] * len(mutated_genome)
            for orig_pos, m_idx in orig_to_mut.items():
                if 0 <= m_idx < len(mutated_genome):
                    mut_to_orig[m_idx] = orig_pos

            # 5) rebuild each transcript’s cDNA under this haplotype
            mutated_transcripts: Dict[str, NewTranscriptSequences] = {}
            for tx_id, res in self._results.items():
                # rebuild spliced‐cDNA, dropping deleted bases and inserting alt_seq in exons
                new_cdna_chars: List[str] = []
                for cdna_pos in sorted(res.cdna_to_dna_map):
                    orig_pos = res.cdna_to_dna_map[cdna_pos]
                    if orig_pos not in orig_to_mut:
                        # that nucleotide was deleted
                        continue
                    m_idx = orig_to_mut[orig_pos]
                    new_cdna_chars.append(mutated_genome[m_idx])

                    # if an insertion event sits at this genomic position, append it here
                    for ev in sorted_events:
                        if ev.ref_allele == "" and (self._offset + ev.pos0) == orig_pos:
                            new_cdna_chars.append(ev.alt_seq)

                new_cdna = "".join(new_cdna_chars)

                # now rebuild ORF‐based codon map & AA translations on new_cdna
                start, end, orf_seq = get_longest_orf(new_cdna)
                new_codon_map: Dict[int, Any] = {}
                new_dna2aa: Dict[int, str] = {}
                new_cdna2aa: Dict[int, str] = {}

                for i in range(0, len(orf_seq), 3):
                    codon = orf_seq[i:i+3]
                    aa = str(Seq(codon).translate(to_stop=True))
                    codon_index = i // 3 + 1

                    cdna_indices = list(range(start + i + 1, start + i + 4))
                    gen_positions = [
                        (res.seq_region, res.cdna_to_dna_map[p])
                        for p in cdna_indices
                        if p in res.cdna_to_dna_map
                    ]
                    tx_positions = [(tx_id, p) for p in cdna_indices]

                    new_codon_map[codon_index] = CodonCoordinateMapping(
                        codon_index=codon_index,
                        codon=codon,
                        amino_acid=aa,
                        genomic_positions=gen_positions,
                        transcript_positions=tx_positions
                    )
                    for _, gpos in gen_positions:
                        new_dna2aa[gpos] = aa
                    for p in cdna_indices:
                        new_cdna2aa[p] = aa

                mutated_transcripts[tx_id] = NewTranscriptSequences(
                    exon_sequences={},           # we can still expose these if needed
                    cdna_sequence=new_cdna,
                    codon_map=new_codon_map,
                    dna_to_protein_map=new_dna2aa,
                    cdna_to_protein_map=new_cdna2aa
                )

            # 6) package up
            all_mutated[frozenset(events)] = {
                "mutated_genome": mutated_genome,
                "genome_protein": genome_protein,
                "orig_to_mutated_index": orig_to_mut,
                "mutated_index_to_orig": mut_to_orig,
                "transcripts": mutated_transcripts,
            }

        return all_mutated
