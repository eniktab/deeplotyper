import logging
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Optional, Sequence
from .data_models import (
    HaplotypeEvent
)
import pysam

logger = logging.getLogger(__name__)


# Type aliases for readability
Pattern = tuple[HaplotypeEvent, ...]
PatternToSamples = dict[Pattern, tuple[str, ...]]


class HaplotypeGroups:
    """Groups of samples sharing identical haplotype event patterns."""

    __slots__ = ("_pattern_to_samples", "_ref_seq", "_start")

    def __init__(
        self,
        pattern_to_samples: PatternToSamples,
        ref_seq: str,
        start: int,
    ) -> None:
        self._pattern_to_samples = pattern_to_samples
        self._ref_seq = ref_seq
        self._start = start

    @property
    def ref_seq(self) -> str:
        """Reference sequence for this haplotype window."""
        return self._ref_seq

    @property
    def start(self) -> int:
        """Zero-based start position of the reference window."""
        return self._start

    @lru_cache(maxsize=None)
    def build_sequence(self, pattern: Pattern) -> str:
        """Build the full haplotype sequence by splicing in events.

        Each event is clamped to the reference bounds and its alternate
        sequence inserted (handling insertions and deletions).
        """
        seq_parts: list[str] = []
        cursor = 0
        ref_len = len(self._ref_seq)

        for ev in sorted(pattern, key=lambda e: e.pos0):
            ev_start = max(0, ev.pos0)
            ev_end = min(ref_len, ev.pos0 + len(ev.ref_allele))

            seq_parts.append(self._ref_seq[cursor:ev_start])
            seq_parts.append(ev.alt_seq)
            cursor = ev_end

        seq_parts.append(self._ref_seq[cursor:])
        return "".join(seq_parts)

    def materialize(self) -> dict[str, tuple[str, ...]]:
        """Return mapping from haplotype sequence to sample tuples."""
        return {
            self.build_sequence(pattern): samples
            for pattern, samples in self._pattern_to_samples.items()
        }

    @classmethod
    def from_vcf(
        cls,
        vcf_path: Path | str,
        ref_seq: str,
        contig: str,
        start: int,
        end: Optional[int] = None,
        samples: Sequence[str] | None = None,
    ) -> "HaplotypeGroups":
        """Read a VCF and create HaplotypeGroups for the given window.

        Args:
            vcf_path (Path | str): Path to the VCF file.
            ref_seq (str): Reference sequence for the window.
            contig (str): Contig name in the VCF.
            start (int): 0-based start position of the window.
            end (Optional[int]): 0-based inclusive end position;
                defaults to start + len(ref_seq) - 1.
            samples (Sequence[str] | None): Subset of sample names to include.

        Raises:
            ValueError: If specified samples are missing from the VCF header.

        Returns:
            HaplotypeGroups: Instance grouping samples by haplotype patterns.
        """
        if end is None:
            end = start + len(ref_seq) - 1

        with pysam.VariantFile(vcf_path) as vcf:
            all_samples = list(vcf.header.samples)
            if samples is None:
                samples = all_samples
            else:
                missing = set(samples) - set(all_samples)
                if missing:
                    raise ValueError(f"Samples not in VCF header: {missing!r}")

            var_patterns: dict[tuple[str, int], list[HaplotypeEvent]] = {
                (sample, hap): [] for sample in samples for hap in (0, 1)
            }

            for rec in vcf.fetch(contig, start, end):
                pos0 = rec.pos - start
                if pos0 < 0 or pos0 >= len(ref_seq):
                    logger.debug(
                        "Skipping variant at %s:%d (pos0=%d) outside window",
                        contig, rec.pos, pos0
                    )
                    continue

                alts = rec.alts or ()
                for sample in samples:
                    call = rec.samples[sample]
                    gt = call.get("GT")
                    if not gt or None in gt:
                        logger.debug(
                            "Skipping missing genotype for sample %s at %s:%d",
                            sample, contig, rec.pos
                        )
                        continue

                    for hap_idx, allele in enumerate(gt):
                        if allele <= 0:
                            continue

                        if allele - 1 >= len(alts):
                            logger.debug(
                                "Allele index %d out of range for record %s:%d",
                                allele, contig, rec.pos
                            )
                            continue

                        alt_seq = alts[allele - 1]
                        if alt_seq.startswith("<") and alt_seq.endswith(">"):
                            logger.debug(
                                "Skipping symbolic allele %s at %s:%d",
                                alt_seq, contig, rec.pos
                            )
                            continue

                        ev = HaplotypeEvent(pos0, rec.ref, alt_seq)
                        var_patterns[(sample, hap_idx)].append(ev)

        tmp: dict[Pattern, set[str]] = defaultdict(set)
        for (sample, _), events in var_patterns.items():
            key = tuple(sorted(events, key=lambda e: e.pos0))
            tmp[key].add(sample)

        pattern_to_samples: PatternToSamples = {
            pat: tuple(sorted(samps))
            for pat, samps in tmp.items()
        }

        return cls(pattern_to_samples, ref_seq, start)