# deeplotyper

[![CI & Release](https://github.com/eniktab/deeplotyper/actions/workflows/ci.yml/badge.svg)](https://github.com/eniktab/deeplotyper/actions/workflows/ci.yml)
[![Publish Python Package](https://github.com/eniktab/deeplotyper/actions/workflows/python-publish.yml/badge.svg)](https://github.com/eniktab/deeplotyper/actions/workflows/python-publish.yml)
[![PyPI version](https://img.shields.io/pypi/v/deeplotyper.svg)](https://pypi.org/project/deeplotyper/)
[![Documentation Status](https://readthedocs.org/projects/deeplotyper/badge/?version=latest)](https://deeplotyper.readthedocs.io/en/latest/)
[//]: # ([![PyPI version](https://img.shields.io/pypi/v/deeplotyper.svg)](https://pypi.org/project/deeplotyper/))

# deeplotyper

Tools for mapping and remapping genomic ↔ transcript sequences.

Deeplotyper is a Python toolkit for genomic and transcriptomic sequence analysis that focuses on mapping coordinates between genomes and transcripts, applying variant haplotypes (sets of SNVs/indels) to reference sequences, and extracting open reading frames (ORFs) as either linear sequences or graph representations. It is designed as an academically rigorous, transparent alternative to traditional variant effect prediction tools. Deeplotyper’s core modules enable fine-grained control and interpretation of complex genetic variants without reliance on large external databases or opaque heuristics.

## Installation

```bash
pip install deeplotyper
```

**Requires:**

- Python ≥ 3.8
- Biopython
- pysam

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
    genome_metadata={"seq_region_accession": "chr1", "start": 100, "strand": 1},
    full_genomic_sequence="ATGGGGTTTCCC...",
    exon_definitions_by_transcript={
        "tx1": [
            {"exon_number": 1, "start": 100, "end": 102, "sequence": "ATG"},
            …
        ]
    },
    transcript_sequences={"tx1": "ATGCCC"},
    exon_orders={"tx1": [1]},
    min_block_length=5
)

# 2. Apply SNV/indel haplotypes
hap_map = {
    (
        HaplotypeEvent(pos0=2, ref_allele="A", alt_seq="G"),
    ): ()
}
remapper = HaplotypeRemapper("ATGAAA...", results)
mutated = remapper.apply_haplotypes(hap_map)

# 3. Group samples by haplotype from a VCF
groups = HaplotypeGroups.from_vcf(
    "variants.vcf.gz",
    ref_seq="ATGAAA...",
    contig="1",
    start=0
)
distinct = groups.materialize()
```

## Sequence Coordinate Mapping (SequenceCoordinateMapper)

One foundational feature of Deeplotyper is coordinate mapping between genomic DNA and transcript (cDNA/mRNA) coordinates. The **SequenceCoordinateMapper** class constructs an internal mapping between a reference sequence (e.g. a genomic region) and one or more transcript definitions (exons/introns structure). This allows conversion of coordinates in both directions (genome → transcript and transcript → genome).

For example, given a gene’s reference DNA sequence and exon coordinates for multiple transcripts (splice variants), the mapper can:

- Translate a genomic position to a position within a transcript (cDNA coordinate).
- Identify which exon or intron a mutation falls into.
- Account for strand orientation and splicing (including reverse-complement mappings).

By building a precise base-level map of exonic regions, SequenceCoordinateMapper provides the groundwork for consistent variant placement across transcripts and enables downstream analyses like coding sequence extraction.

**Implementation detail:** Internally, the mapper may produce a linear coordinate index for each transcript relative to the reference. For instance, if Transcript A has exons 1–100 and 201–300 on the reference genome, a coordinate like genomic 250 can be mapped to position 150 of Transcript A’s cDNA.

## Haplotype Remapping (HaplotypeRemapper)

Deeplotyper supports applying a set of genetic variants — collectively forming a haplotype — onto reference sequences or transcripts. The **HaplotypeRemapper** class takes a SequenceCoordinateMapper and a haplotype map (a collection of variants such as SNVs, insertions, deletions, or complex multi-nucleotide changes) and remaps the reference sequence to produce the altered (haplotype) sequence.

- Ensures variants are applied in the correct positions across multi-exon transcripts.
- Handles insertions and deletions (indels), adjusting downstream coordinates.
- Supports complex events like multi-base substitutions or combinations of proximal variants.
- Can model gene fusions or structural rearrangements by mapping coordinates from two reference sequences into one combined transcript.

The output is typically a new sequence (e.g. the mutated cDNA), along with diffs or lists of changed positions for full transparency.

## ORF Extraction (find_orfs and get_longest_orf)

To assess coding impacts, Deeplotyper can extract open reading frames (ORFs) from sequences:

- **find_orfs** scans a nucleotide sequence to identify all ORFs bounded by start and stop codons in the correct reading frame.
- **get_longest_orf** retrieves the longest ORF from a given sequence.

These functions help reveal variant-induced effects such as novel start codons, truncated proteins, or frameshifts. Graph representations of ORFs (nodes = exons/segments, edges = splice connections) are also supported for visualizing complex haplotypes.

## Sequence Alignment (make_aligner and apply_alignment_gaps)

When visualizing indels, Deeplotyper provides utilities for pairwise sequence alignment:

- **make_aligner** returns a configured Biopython PairwiseAligner (global or local modes).
- **apply_alignment_gaps** projects alignment gaps onto coordinate mappings or sequence strings, inserting dashes (‐) to show indels.

Example alignment output:

```
Ref: ATGCCCACGT...
Alt: ATG--ACGT...
```

This aids in interpreting frameshifts or in-frame indels and their effects on codon numbering.

## Linear Coordinate Construction (build_linear_coords)

The **build_linear_coords** utility flattens a spliced transcript into a continuous cDNA or protein coordinate space and maps it back to genomic coordinates. Useful for:

- Creating lookup tables (e.g. transcript→genome).
- Plotting gene models.
- Adjusting coordinates after indels in haplotype transcripts.

## Example Use Case

```python
from deeplotyper import SequenceCoordinateMapper, HaplotypeRemapper, find_orfs, get_longest_orf

# 1. Reference sequence (toy example)
gene_name = "GENE1"
chrom = "chr1"
strand = "+"

reference_seq = (
    "ATGGTcacct...TTAG"
)

# 2. Exon definitions for two transcripts
transcript1_exons = [(1, 300), (401, 600)]
transcript2_exons = [(1, 300), (501, 700)]

transcripts = {
    "Transcript1": {"exons": transcript1_exons, "strand": "+", "cds_start": 1, "cds_end": 600},
    "Transcript2": {"exons": transcript2_exons, "strand": "+", "cds_start": 1, "cds_end": 700}
}

mapper = SequenceCoordinateMapper(reference_seq, transcripts)

# 3. Define a haplotype (list of variant dicts)
haplotype = [
    {"pos": 50,  "ref": "G",   "alt": "T"},
    {"pos": 310, "ref": "",    "alt": "ACG"},
    {"pos": 450, "ref": "AGCT","alt": ""},
    {"pos": 480, "ref": "A",   "alt": "TT"},
]

remapper = HaplotypeRemapper(mapper, haplotype)

mut_seq_t1 = remapper.get_sequence("Transcript1")
mut_seq_t2 = remapper.get_sequence("Transcript2")

print(f"Transcript1 (mutated) length: {len(mut_seq_t1)}")
print(mut_seq_t1[40:60])

# 5. ORF extraction in mutated Transcript1
orfs = find_orfs(mut_seq_t1, assume_start_codon=True)
longest_orf = get_longest_orf(mut_seq_t1)
print(f"Number of ORFs: {len(orfs)}")
print(f"Longest ORF length: {len(longest_orf)}")
```

## Addressing Limitations of VEP and Haplosaurus

Traditional VEP/Haplosaurus workflows have known limitations:

- **Complex variant support:** Doesn’t natively handle gene fusions, multi-exon deletions, or intronic/splice-site changes. Deeplotyper applies any user-specified set of variants.
- **Database dependency:** Requires multi-GB Ensembl caches and compiled APIs. Deeplotyper is pure-Python and works on user-provided sequences/coords.
- **Edge cases:** Can fail on short transcripts or produce opaque “high impact” labels. Deeplotyper’s transparent implementation traces frameshifts and disrupted sequences.
- **Opacity:** VEP uses black-box predictors (SIFT/PolyPhen). Deeplotyper exposes explicit sequence changes, enabling direct inspection of altered codons or ORFs.

## License

[MIT](LICENSE)

## Contributing
We welcome contributions! Feel free to add requests in the issues section or directly contribute with a pull request.

## Citations

- [Haplosaurus computes protein haplotypes for use in precision drug design | Nature Communications](https://www.nature.com/articles/s41467-018-06542-1)
- [GitHub - Ensembl/ensembl-vep: The Ensembl Variant Effect Predictor predicts the functional effects of genomic variants](https://github.com/Ensembl/ensembl-vep)
- [Haplosaurus can require more than 60GB memory for a single ...](https://github.com/Ensembl/ensembl-vep/issues/497)
- [Gene-specific features enhance interpretation of mutational impact on acid alpha-glucosidase enzyme activity - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7329270/)
- [The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity - PubMed](https://pubmed.ncbi.nlm.nih.gov/25684150/)