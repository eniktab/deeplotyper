from typing import Dict, List, Optional, Tuple
from collections import defaultdict
from .data_models import (
    RawBase,
    BaseCoordinateMapping,
    CodonCoordinateMapping,
    SequenceMappingResult,
    TranscriptMappingResult,
    ExonDef,
    MappingEntry,
)
from .alignment_utils import make_aligner, apply_alignment_gaps
from .coordinate_utils import build_raw_genome_coords, build_raw_transcript_coords
from Bio.Seq import Seq


class SequenceCoordinateMapper:
    """Maps between spliced (cDNA), transcript, and genomic sequences."""

    def _find_exon_alignment_blocks(
        self,
        spliced: str,
        transcript: str,
        exon_order: Optional[List[int]],
        exon_defs: Optional[Dict[int, ExonDef]],
    ) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
        """Return aligned blocks between spliced and transcript sequences.

        Aligns either per-exon (using provided exon_order and exon_defs) or the
        entire spliced vs. transcript sequence (global alignment) to compute
        matching intervals ("blocks") in both strings.

        Args:
            spliced: The concatenated spliced (cDNA) sequence.
            transcript: The full transcript sequence.
            exon_order: Order of exon numbers to align individually.
            exon_defs: Definitions for each exon keyed by exon number.

        Returns:
            A pair of lists:
              - blocks_spliced: list of (start, end) intervals in `spliced`
              - blocks_transcript: corresponding intervals in `transcript`
        """
        local_aln = make_aligner("local")
        blocks_spliced: List[Tuple[int, int]] = []
        blocks_transcript: List[Tuple[int, int]] = []
        cursor_sp, cursor_tx = 0, 0

        if exon_order and exon_defs:
            for ex in exon_order:
                seq = exon_defs[ex]["sequence"]
                length = len(seq)
                pos = transcript.find(seq, cursor_tx)

                if pos != -1:
                    blk_sp = [(cursor_sp, cursor_sp + length)]
                    blk_tx = [(pos, pos + length)]
                else:
                    aln = local_aln.align(seq, transcript)[0]
                    raw_sp, raw_tx = aln.aligned
                    blk_sp = [(cursor_sp + int(a), cursor_sp + int(b)) for a, b in raw_sp]
                    blk_tx = [(int(a), int(b)) for a, b in raw_tx]

                blocks_spliced.extend(blk_sp)
                blocks_transcript.extend(blk_tx)
                cursor_sp += length
                cursor_tx = max(end for _, end in blk_tx)
        else:
            global_aln = make_aligner("global")
            aln = global_aln.align(spliced, transcript)[0]
            raw_sp, raw_tx = aln.aligned
            # Convert to Python ints and tuples to avoid numpy types
            blocks_spliced = [(int(a), int(b)) for a, b in raw_sp]
            blocks_transcript = [(int(a), int(b)) for a, b in raw_tx]

        return blocks_spliced, blocks_transcript

    def _extract_base_mappings(
        self,
        spliced: str,
        transcript: str,
        raw_genome: List[RawBase],
        raw_transcript: List[RawBase],
        blocks_sp: List[Tuple[int, int]],
        blocks_tx: List[Tuple[int, int]],
        min_block_length: int
    ) -> List[BaseCoordinateMapping]:
        """Flatten aligned blocks into BaseCoordinateMapping, filtering out short runs.

        Applies gap insertion on the spliced vs. transcript strings, finds
        contiguous runs of base matches, filters runs shorter than
        `min_block_length`, and converts each matched base pair into a
        BaseCoordinateMapping.

        Args:
            spliced: The spliced (cDNA) sequence.
            transcript: The transcript sequence.
            raw_genome: List of RawBase for the genome coordinates.
            raw_transcript: List of RawBase for the transcript coordinates.
            blocks_sp: Aligned block intervals in `spliced`.
            blocks_tx: Aligned block intervals in `transcript`.
            min_block_length: Minimum length of a contiguous match run to include.

        Returns:
            A list of BaseCoordinateMapping entries for each retained base match.
        """
        sp_g, tx_g = apply_alignment_gaps(spliced, transcript, blocks_sp, blocks_tx)

        # collect matches
        matches: List[Tuple[int, int, RawBase, RawBase]] = []
        i_sp = i_tx = 0
        for a_char, t_char in zip(sp_g, tx_g):
            if a_char != "-" and t_char != "-":
                matches.append((i_sp, i_tx, raw_genome[i_sp], raw_transcript[i_tx]))
            if a_char != "-":
                i_sp += 1
            if t_char != "-":
                i_tx += 1

        # group into consecutive runs
        runs: List[List[Tuple[int, int, RawBase, RawBase]]] = []
        if matches:
            current = [matches[0]]
            for sp_i, tx_i, g_raw, t_raw in matches[1:]:
                _, _, prev_g, prev_t = current[-1]
                if (g_raw.position == prev_g.position + 1
                        and t_raw.position == prev_t.position + 1):
                    current.append((sp_i, tx_i, g_raw, t_raw))
                else:
                    runs.append(current)
                    current = [(sp_i, tx_i, g_raw, t_raw)]
            runs.append(current)

        # emit mappings for runs ≥ min_block_length
        result: List[BaseCoordinateMapping] = []
        for run in runs:
            if len(run) < min_block_length:
                continue
            for _, _, g_raw, t_raw in run:
                result.append(BaseCoordinateMapping(
                    global_position     = g_raw.position,
                    genomic_position    = g_raw.position,
                    transcript_position = t_raw.position,
                    genomic_base        = g_raw.base,
                    transcript_base     = t_raw.base
                ))
        return result

    def _generate_base_coordinate_mappings(
        self,
        spliced_sequence: str,
        genome_metadata: Dict[str, any],
        exon_order: Optional[List[int]],
        exon_definitions: Optional[List[ExonDef]],
        transcript_sequence: str,
        full_genomic_sequence: str,
        transcript_start: int,
        transcript_name: str,
        min_block_length: int
    ) -> SequenceMappingResult:
        """Produce base-level mappings and gapped sequences per exon.

        Builds raw coordinate lists for genome and transcript, finds aligned
        blocks, performs per-exon full-genome to exon alignments, and extracts
        base mappings.

        Args:
            spliced_sequence: Concatenated exon sequences (cDNA).
            genome_metadata: Metadata including 'seq_region_accession',
                optionally 'strand' and 'start'.
            exon_order: Ordered list of exon numbers.
            exon_definitions: Exon definitions list.
            transcript_sequence: Full transcript sequence.
            full_genomic_sequence: Entire genomic region sequence.
            transcript_start: Starting index for transcript positions.
            transcript_name: Identifier for the transcript.
            min_block_length: Minimum contiguous match length for base mapping.

        Returns:
            SequenceMappingResult containing base mappings, sorted mappings,
            gapped sequences, and alignment maps.
        """
        region = genome_metadata["seq_region_accession"]
        strand = genome_metadata.get("strand", 1)
        offset = genome_metadata.get("start", 1)
        exon_defs = {e["exon_number"]: e for e in (exon_definitions or [])}

        blocks_sp, blocks_tx = self._find_exon_alignment_blocks(
            spliced_sequence, transcript_sequence,
            exon_order, exon_defs
        )

        raw_genome = build_raw_genome_coords(
            spliced_sequence, region, strand, offset, exon_order, exon_defs
        )
        raw_transcript = build_raw_transcript_coords(
            transcript_sequence, transcript_name, transcript_start
        )

        # full-genome-per-exon alignment
        g_full, g_exon, g_map = {}, {}, {}
        for exnum in (exon_order or sorted(exon_defs)):
            seq = exon_defs[exnum]["sequence"]
            aln = make_aligner("local").align(full_genomic_sequence, seq)[0]
            b_g, b_e = aln.aligned
            gapped_full, gapped_exon = apply_alignment_gaps(
                full_genomic_sequence, seq, b_g, b_e
            )
            g_full[exnum] = gapped_full
            g_exon[exnum] = gapped_exon

            pairs: List[Tuple[int, int]] = []
            idx_g = idx_e = 0
            for cg, ce in zip(gapped_full, gapped_exon):
                if cg != "-" and ce != "-":
                    pairs.append((offset + idx_g, idx_e + 1))
                if cg != "-":
                    idx_g += 1
                if ce != "-":
                    idx_e += 1
            g_map[exnum] = pairs

        base_maps = self._extract_base_mappings(
            spliced_sequence,
            transcript_sequence,
            raw_genome,
            raw_transcript,
            blocks_sp, blocks_tx,
            min_block_length
        )
        sorted_by_genome = sorted(base_maps, key=lambda m: m.genomic_position)

        return SequenceMappingResult(
            base_mappings                   = base_maps,
            base_mappings_by_genome         = sorted_by_genome,
            gapped_full_genome_sequence     = g_full,
            gapped_transcript_sequence      = g_exon,
            gapped_alignment_map            = g_map
        )

    def _generate_codon_mappings(
        self,
        base_maps: List[BaseCoordinateMapping],
        seq_region: str,
        transcript_id: str
    ) -> Tuple[Dict[int, CodonCoordinateMapping], Dict[int, str], Dict[int, str]]:
        """Build codon to amino-acid mappings from base-level mappings.

        Groups every three consecutive BaseCoordinateMapping entries into a
        codon, translates to amino acid, and constructs lookup maps.

        Args:
            base_maps: List of base mappings.
            seq_region: Sequence region accession for genomic positions.
            transcript_id: Identifier for the transcript.

        Returns:
            A tuple of:
              - codon_map: maps codon index to CodonCoordinateMapping
              - dna2aa: maps genomic position to amino acid
              - cdna2aa: maps transcript position to amino acid
        """
        raw_g = [(seq_region, m.genomic_position, m.genomic_base) for m in base_maps]
        raw_t = [(transcript_id, m.transcript_position, m.transcript_base) for m in base_maps]

        codon_map: Dict[int, CodonCoordinateMapping] = {}
        dna2aa: Dict[int, str] = {}
        cdna2aa: Dict[int, str] = {}
        for i in range(0, len(raw_g) - 2, 3):
            trip_g = raw_g[i : i + 3]
            trip_t = raw_t[i : i + 3]
            cod = "".join(b for _, _, b in trip_g)
            aa = str(Seq(cod).translate(to_stop=False))
            idx = i // 3 + 1

            codon_map[idx] = CodonCoordinateMapping(
                codon_index          = idx,
                codon                = cod,
                amino_acid           = aa,
                genomic_positions    = [(r, p) for r, p, _ in trip_g],
                transcript_positions = [(r, p) for r, p, _ in trip_t]
            )
            for _, p, _ in trip_g:
                dna2aa[p] = aa
            for _, p, _ in trip_t:
                cdna2aa[p] = aa

        return codon_map, dna2aa, cdna2aa

    def map_transcripts(
        self,
        genome_metadata: Dict[str, any],
        full_genomic_sequence: str,
        exon_definitions_by_transcript: Dict[str, List[ExonDef]],
        transcript_sequences: Dict[str, str],
        exon_orders: Dict[str, List[int]],
        min_block_length: int = 15,
    ) -> Dict[str, TranscriptMappingResult]:
        """Public API: map each transcript to genomic coordinates, codons, and amino acids.

        Iterates over transcripts, constructs spliced sequences from exon definitions,
        generates base and codon mappings, and assembles a TranscriptMappingResult.

        Args:
            genome_metadata: Metadata with keys 'seq_region_accession',
                optionally 'strand' and 'start'.
            full_genomic_sequence: Entire genomic region sequence.
            exon_definitions_by_transcript: ExonDef lists per transcript.
            transcript_sequences: Raw transcript sequences by transcript ID.
            exon_orders: Exon order lists per transcript ID.
            min_block_length: Minimum contiguous match length for base mapping.

        Returns:
            A dict mapping each transcript ID to its TranscriptMappingResult.
        """
        meta = genome_metadata.copy()
        meta.setdefault("seq_region_accession", meta.get("seq_region_name"))
        results: Dict[str, TranscriptMappingResult] = {}

        for tx_id, tx_seq in transcript_sequences.items():
            defs = exon_definitions_by_transcript.get(tx_id)
            if not defs:
                continue
            order = exon_orders.get(tx_id, sorted(e["exon_number"] for e in defs))
            by_num = {e["exon_number"]: e for e in defs}
            spliced = "".join(by_num[n]["sequence"] for n in order)

            seq_res = self._generate_base_coordinate_mappings(
                spliced_sequence       = spliced,
                genome_metadata        = meta,
                exon_order             = order,
                exon_definitions       = defs,
                transcript_sequence    = tx_seq,
                full_genomic_sequence  = full_genomic_sequence,
                transcript_start       = 1,
                transcript_name        = tx_id,
                min_block_length       = min_block_length
            )

            base_map: List[MappingEntry] = [
                MappingEntry(
                    transcript_index = m.transcript_position,
                    transcript_base  = m.transcript_base,
                    genomic_index    = m.genomic_position,
                    genomic_base     = m.genomic_base
                )
                for m in seq_res.base_mappings
            ]

            cdna2dna = {e["transcript_index"]: e["genomic_index"] for e in base_map}
            dna2cdna = {e["genomic_index"]: e["transcript_index"] for e in base_map}

            exon_to_genome: Dict[int, List[MappingEntry]] = defaultdict(list)
            exon_to_tx: Dict[int, List[MappingEntry]] = defaultdict(list)
            dna_to_exon: Dict[int, int] = {}
            cdna_to_exon: Dict[int, int] = {}
            cur = 1
            for ex in order:
                seq = by_num[ex]["sequence"]
                for i in range(len(seq)):
                    cdna_to_exon[cur + i] = ex
                cur += len(seq)

            for e in base_map:
                pos_g, pos_t = e["genomic_index"], e["transcript_index"]
                ex = cdna_to_exon.get(pos_t)
                if ex is not None:
                    exon_to_tx[ex].append(e)
                if any(pos_g == p for p in range(by_num[ex]["start"], by_num[ex]["end"] + 1)):
                    exon_to_genome[ex].append(e)

            codon_map, dna2aa, cdna2aa = self._generate_codon_mappings(
                seq_res.base_mappings,
                meta["seq_region_accession"],
                tx_id
            )

            results[tx_id] = TranscriptMappingResult(
                transcript_to_genomic       = base_map,
                cdna_to_dna_map             = cdna2dna,
                dna_to_cdna_map             = dna2cdna,
                exon_to_genomic             = dict(exon_to_genome),
                exon_to_transcript          = dict(exon_to_tx),
                dna_to_exon_map             = dna_to_exon,
                cdna_to_exon_map            = cdna_to_exon,
                codon_map                   = codon_map,
                dna_to_protein_map          = dna2aa,
                cdna_to_protein_map         = cdna2aa,
                gapped_full_genome_sequence = seq_res.gapped_full_genome_sequence,
                gapped_transcript_sequence  = seq_res.gapped_transcript_sequence,
                gapped_alignment_map        = seq_res.gapped_alignment_map,
                seq_region                  = meta["seq_region_accession"],
                offset                      = meta.get("start", 1),
            )

        return results
