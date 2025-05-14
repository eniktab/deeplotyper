from typing import (
    Tuple,List
)

def find_orfs(seq: str) -> List[Tuple[int, int, str]]:
    """
    Identify all open reading frames (ORFs) in the given nucleotide sequence.
    Searches only the forward strand in reading frames 0, 1, and 2.
    """
    seq = seq.upper()
    orfs: List[Tuple[int,int,str]] = []
    for frame in range(3):
        start_index = None
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if start_index is None:
                if codon == "ATG":
                    start_index = i
            else:
                if codon in ("TAA", "TAG", "TGA"):
                    orfs.append((start_index, i+3, seq[start_index:i+3]))
                    start_index = None
    return orfs

def get_longest_orf(seq: str) -> Tuple[int, int, str]:
    """
    Return the longest ORF found in the sequence, as (start, end, orf_seq).
    If no ORF is found, treat the entire sequence as the ORF.
    """
    orfs = find_orfs(seq)
    if not orfs:
        return (0, len(seq), seq)
    return max(orfs, key=lambda x: len(x[2]))