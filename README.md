# deeplotyper

[![CI & Release](https://github.com/eniktab/deeplotyper/actions/workflows/python-publish.yml/badge.svg)](https://github.com/eniktab/deeplotyper/actions/workflows/ci.yml)
[![PyPI version](https://img.shields.io/pypi/v/deeplotyper.svg)](https://pypi.org/project/deeplotyper/)


Tools for graph based and string based mapping and remapping genomic ↔ transcript ↔ aminoacid sequences.

## Installation

```bash
pip install deeplotyper
```

> **Requires**  
> - Python ≥ 3.8  
> - Biopython  
> - pysam  

## Quickstart

```python
from deeplotyper import (
    SequenceCoordinateMapper,
    HaplotypeRemapper,
    HaplotypeGroups,
    find_orfs, get_longest_orf,
    make_aligner, apply_alignment_gaps,
    build_linear_coords, build_raw_genome_coords, build_raw_transcript_coords,
    BaseCoordinateMapping, CodonCoordinateMapping,
    SequenceMappingResult, TranscriptMappingResult,
    HaplotypeEvent, NewTranscriptSequences, RawBase
)

# 1. Map a transcript to the genome
mapper = SequenceCoordinateMapper()
results = mapper.map_transcripts(
    genome_metadata={ "seq_region_accession": "chr1", "start": 100, "strand": 1 },
    full_genomic_sequence="ATGGGGTTTCCC...",
    exon_definitions_by_transcript={ "tx1": [ { "exon_number":1, "start":100, "end":102, "sequence":"ATG" }, … ] },
    transcript_sequences={ "tx1":"ATGCCC" },
    exon_orders={ "tx1":[1] },
    min_block_length=5
)

# 2. Apply SNV/indel haplotypes
hap_map = {
    ( HaplotypeEvent(pos0=2, ref_allele="A", alt_seq="G"), ): ()
}
remapper = HaplotypeRemapper("ATGAAA...", results)
mutated = remapper.apply_haplotypes(hap_map)

# 3. Group samples by haplotype from a VCF
groups = HaplotypeGroups.from_vcf("variants.vcf.gz", ref_seq="ATGAAA...", contig="1", start=0)
distinct = groups.materialize()
```  

## Package Modules

### `orf_utils`  
Functions to detect and extract open reading frames.  
- `find_orfs(seq: str) → List[(start,end,seq)]`  
- `get_longest_orf(seq: str) → (start,end,orf_seq)`

### `alignment_utils`  
Wrappers around Biopython’s `PairwiseAligner` for caching and gap insertion.  
- `make_aligner(...) → PairwiseAligner`  
- `apply_alignment_gaps(seq1, seq2, blocks1, blocks2) → (gapped1, gapped2)`

### `data_models`  
Typed dataclasses for coordinate‐mapping results and events.  
- `BaseCoordinateMapping`, `CodonCoordinateMapping`  
- `SequenceMappingResult`, `TranscriptMappingResult`  
- `RawBase`, `HaplotypeEvent`, `NewTranscriptSequences`

### `coordinate_utils`  
Builds raw coordinate lists from sequences or exon definitions.  
- `build_linear_coords(sequence, region, start_pos, strand)`  
- `build_raw_genome_coords(spliced, region, strand, start_offset, exon_order=None, exon_defs=None)`  
- `build_raw_transcript_coords(transcript, region, start_pos)`

### `mapper`  
`SequenceCoordinateMapper` — maps exons/transcripts to genomic coordinates at base and codon resolution, producing both gapped alignments and mapping objects.

### `remapper`  
`HaplotypeRemapper` — given a reference slice and transcript mappings, applies SNV/indel haplotypes to rebuild mutated genome & transcript sequences, ORFs, and all coordinate maps.

### `vcf_haploevents`  
`HaplotypeGroups` — read a VCF window to group samples by identical haplotype event patterns, and materialize haplotype→sample mappings. 

## Contributing

1. Fork the repo  
2. Create a feature branch (`git checkout -b feat/…`)  
3. Add tests under `tests/` (pytest)  
4. Submit a pull request  

## License

[MIT](LICENSE)
