import logging
from functools import lru_cache
from typing import Any, Dict, List, Optional, Sequence, Tuple

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


class NCBIGeneSequenceFetcher:
    """Fetch gene‐centric sequences from NCBI Entrez.

    Provides methods to retrieve genomic regions, cDNA, protein sequences,
    and exon information for a specified gene in a given organism and
    assembly, with caching to avoid redundant calls.
    """

    _ASSEMBLY_MAPPING: Dict[str, str] = {
        "GRCh38.p14": "GCF_000001405.40",
        "T2T-CHM13v2.0": "GCA_009914755.3",
    }

    def __init__(
        self,
        email: Optional[str] = None,
        organism: str = "Homo sapiens",
        assembly: str = "GRCh38.p14",
    ) -> None:
        """Initialize fetcher with Entrez email, organism, and assembly.

        Args:
            email: Email for Entrez; uses placeholder if None.
            organism: Organism name for gene queries.
            assembly: Assembly key; must be in _ASSEMBLY_MAPPING.

        Raises:
            ValueError: If assembly is not supported.
        """
        Entrez.email = email or "anonymous@example.com"
        if email is None:
            logger.warning(
                "Entrez email not provided; using placeholder address."
            )

        if assembly not in self._ASSEMBLY_MAPPING:
            supported = ", ".join(self._ASSEMBLY_MAPPING)
            raise ValueError(
                f"Assembly '{assembly}' is not supported. "
                f"Choose one of: {supported}"
            )

        self.organism = organism
        self.assembly = assembly
        self._assembly_acc = self._ASSEMBLY_MAPPING[assembly]

    def fetch_gene_id(self, gene_name: str) -> str:
        """Retrieve the NCBI Gene ID for a given gene name.

        Args:
            gene_name: Gene symbol or full name.

        Returns:
            The first matching NCBI Gene ID.

        Raises:
            ValueError: If no matching gene is found.
        """
        term = (
            f'{gene_name}[Gene Name] AND "{self.organism}"[Organism]'
        )
        with Entrez.esearch(db="gene", term=term, retmax=1) as handle:
            record = Entrez.read(handle)

        ids = record.get("IdList", [])
        if not ids:
            raise ValueError(
                f"No gene found for '{gene_name}' in '{self.organism}'."
            )
        return ids[0]

    def fetch_gene_summary(self, gene_id: str) -> Dict[str, Any]:
        """Fetch summary information for a gene ID.

        Args:
            gene_id: NCBI Gene ID.

        Returns:
            Parsed Entrez summary as a dictionary.
        """
        with Entrez.esummary(db="gene", id=gene_id, retmode="xml") as handle:
            return Entrez.read(handle)

    def fetch_genomic_sequence(
        self, accession: str, start: int, stop: int
    ) -> SeqRecord:
        """Fetch a genomic sequence region by accession and coordinates.

        Args:
            accession: Accession of the genomic record.
            start: Zero-based start (inclusive).
            stop: Zero-based end (inclusive).

        Returns:
            SeqRecord of the requested region.
        """
        with Entrez.efetch(
            db="nuccore",
            id=accession,
            rettype="fasta",
            retmode="text",
            seq_start=start + 1,
            seq_stop=stop + 1,
        ) as handle:
            return SeqIO.read(handle, "fasta")

    def fetch_cdna_sequences(self, gene_id: str) -> List[SeqRecord]:
        """Fetch all RefSeq cDNA transcripts linked to a gene.

        Args:
            gene_id: NCBI Gene ID.

        Returns:
            List of cDNA SeqRecord objects.
        """
        rna_ids = self._fetch_linked_record_ids(
            gene_id, db="nuccore", linkname="gene_nuccore_refseqrna"
        )
        return self._fetch_fasta_records(rna_ids, db="nuccore")

    def fetch_protein_sequences(self, gene_id: str) -> List[SeqRecord]:
        """Fetch all RefSeq protein translations linked to a gene.

        Args:
            gene_id: NCBI Gene ID.

        Returns:
            List of protein SeqRecord objects.
        """
        prot_ids = self._fetch_linked_record_ids(
            gene_id, db="protein", linkname="gene_protein_refseq"
        )
        return self._fetch_fasta_records(prot_ids, db="protein")

    def fetch_all_ncbi_sequences(self, gene_name: str) -> Dict[str, Any]:
        """Fetch comprehensive sequence data for a gene.

        Args:
            gene_name: Gene symbol or full name.

        Returns:
            Dictionary containing genomic, cDNA, MANE, protein, exon data,
            and region information.
        """
        gene_id = self.fetch_gene_id(gene_name)
        summary = self.fetch_gene_summary(gene_id)
        region = self._select_genomic_region_info(summary)

        genomic = self.fetch_genomic_sequence(
            region["ChrAccVer"],
            int(region["ChrStart"]),
            int(region["ChrStop"]),
        )
        cdna = self.fetch_cdna_sequences(gene_id)
        protein = self.fetch_protein_sequences(gene_id)
        exon_info = self.get_exon_info_ncbi(gene_id)
        exon_arrangements = self.get_exon_arrangement_ncbi(gene_id)

        mane_cdna: Optional[SeqRecord] = None
        mane_protein: Optional[SeqRecord] = None
        try:
            mane_cdna = self.fetch_mane_select_cdna(gene_id)
            mane_protein = self.fetch_mane_select_protein(gene_id)
        except ValueError:
            # No annotated or multiple MANE Select transcripts found
            pass

        return {
            "genomic": genomic,
            "cdna": cdna,
            "mane_cdna": mane_cdna,
            "protein": protein,
            "mane_protein": mane_protein,
            "exon_info": exon_info,
            "exon_arrangements": exon_arrangements,
            "region_info": region,
        }

    def fetch_all(self, gene_name: str) -> Dict[str, Any]:
        raw = self.fetch_all_ncbi_sequences(gene_name)

        # Convert SeqRecord objects to plain strings
        cdna_map = {rec.id: str(rec.seq) for rec in raw["cdna"]}
        protein_map = {rec.id: str(rec.seq) for rec in raw["protein"]}

        # region_info → genomic_info, genomic SeqRecord → string
        return {
            "genomic_sequence": str(raw["genomic"].seq),
            "genomic_info": self._normalize_genomic_info(raw["region_info"]),
            "cdna": cdna_map,
            "protein": protein_map,
            # NCBI has no CCDS map
            "ccds": {},
            # pull out MANE Select if present
            "mane_cdna": (str(raw["mane_cdna"].seq)
                          if raw.get("mane_cdna") else None),
            "mane_protein": (str(raw["mane_protein"].seq)
                             if raw.get("mane_protein") else None),
            "exon_info": raw["exon_info"],
            "exon_arrangements": raw["exon_arrangements"],
        }

    def get_exon_info_ncbi(
        self, gene_id: str, transcript_ids: Optional[List[str]] = None
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Retrieve exon coordinates and sequences for transcripts.

        If transcript_ids is None, all RefSeq cDNA transcript IDs are used.

        Args:
            gene_id: NCBI Gene ID.
            transcript_ids: Optional list of transcript IDs to fetch.

        Returns:
            Mapping from transcript ID to list of exon dictionaries.
        """
        if transcript_ids is None:
            transcript_ids = self._fetch_linked_record_ids(
                gene_id, db="nuccore", linkname="gene_nuccore_refseqrna"
            )
        if not transcript_ids:
            return {}

        records = self._fetch_genbank_records(transcript_ids, db="nuccore")
        exon_data: Dict[str, List[Dict[str, Any]]] = {}

        for rec in records:
            exons = [f for f in rec.features if f.type == "exon"]
            if not exons:
                continue

            strand = exons[0].location.strand or 1
            # Sort exons by genomic coordinate, respecting strand
            exons.sort(
                key=lambda f: int(f.location.start),
                reverse=(strand < 0),
            )

            entries: List[Dict[str, Any]] = []
            for idx, feat in enumerate(exons, start=1):
                raw = feat.qualifiers.get("number", [None])[0]
                try:
                    num = int(raw)
                except (TypeError, ValueError):
                    num = idx

                start = int(feat.location.start) + 1
                end = int(feat.location.end)
                seq = str(feat.extract(rec.seq))

                entries.append({
                    "exon_number": num,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "sequence": seq,
                })

            exon_data[rec.id] = entries

        return exon_data

    def get_exon_arrangement_ncbi(
            self, gene_id: str) -> Dict[str, Dict[str, List[int]]]:
        """Determine canonical vs alternative exon arrangements.

        Args:
            gene_id: NCBI Gene ID.

        Returns:
            Dictionary with 'canonical' and 'alternatives' exon-number lists.
        """
        info = self.get_exon_info_ncbi(gene_id)
        return self._compute_exon_arrangement(info)

    def get_mane_select_transcript_id(self, gene_id: str) -> str:
        """Retrieve the MANE Select RefSeq transcript ID for a gene.

        Args:
            gene_id: NCBI Gene ID.

        Returns:
            Single MANE Select transcript ID.

        Raises:
            ValueError: If none or multiple MANE Select transcripts are found.
        """
        tids = self._fetch_linked_record_ids(
            gene_id, db="nuccore", linkname="gene_nuccore_refseqrna"
        )
        recs = self._fetch_genbank_records(tids, db="nuccore")

        mane = {
            rec.id
            for rec in recs
            if any("mane" in kw.lower() for kw in rec.annotations.get("keywords", []))
        }
        if not mane:
            raise ValueError(
                f"No MANE Select transcript found for gene {gene_id}"
            )
        if len(mane) > 1:
            raise ValueError(
                f"Multiple MANE Select transcripts found: {sorted(mane)!r}"
            )
        return mane.pop()

    def fetch_mane_select_cdna(self, gene_id: str) -> SeqRecord:
        """Fetch the full GenBank cDNA record for the MANE Select transcript.

        Args:
            gene_id: NCBI Gene ID.

        Returns:
            SeqRecord for the MANE Select cDNA.
        """
        tid = self.get_mane_select_transcript_id(gene_id)
        records = self._efetch_records(
            (tid,), db="nuccore", rettype="gb", parser="genbank"
        )
        return records[0]

    def fetch_mane_select_protein(self, gene_id: str) -> SeqRecord:
        """Fetch the protein translation for the MANE Select transcript.

        Args:
            gene_id: NCBI Gene ID.

        Returns:
            SeqRecord of the MANE Select protein.

        Raises:
            ValueError: If no CDS/protein is found for the transcript.
        """
        cdna = self.fetch_mane_select_cdna(gene_id)
        for feat in cdna.features:
            if feat.type == "CDS":
                prot_ids = feat.qualifiers.get("protein_id", [])
                if prot_ids:
                    return self._fetch_fasta_records(prot_ids, db="protein")[0]

        raise ValueError(
            f"No CDS/protein found for MANE Select of gene {gene_id}"
        )

    @staticmethod
    def _normalize_genomic_info(raw: Dict[str, Any]) -> Dict[str, Any]:
        """Map NCBI region_info dict into our unified genomic_info schema."""
        # raw keys:
        # 'AnnotationRelease','AssemblyAccVer','ChrAccVer','ChrStart','ChrStop'
        return {
            "assembly": raw.get("AssemblyAccVer"),
            "annotation_release": raw.get("AnnotationRelease"),
            "seq_region_name": None,
            "seq_region_accession": raw.get("ChrAccVer"),
            "start": int(raw.get("ChrStart", 0)),
            "end": int(raw.get("ChrStop", 0))
        }

    def _select_genomic_region_info(
        self, summary: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Pick genomic region info matching target assembly or fallback.

        Args:
            summary: Parsed summary dictionary from Entrez.

        Returns:
            A dict containing region info keys: 'ChrAccVer', 'ChrStart', 'ChrStop'.

        Raises:
            ValueError: If no region information is available.
        """
        doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
        for hist in doc.get("LocationHist", []):
            if hist.get("AssemblyAccVer") == self._assembly_acc:
                return hist

        gen_list = doc.get("GenomicInfo", [])
        if gen_list:
            first = gen_list[0]
            logger.warning(
                "Assembly '%s' not found in LocationHist; using "
                "GenomicInfo accession '%s'.",
                self.assembly,
                first.get("ChrAccVer"),
            )
            return first

        raise ValueError("No genomic region information in summary.")

    @staticmethod
    @lru_cache(maxsize=None)
    def _efetch_records(
        ids: Tuple[str, ...], db: str, rettype: str, parser: str
    ) -> Tuple[SeqRecord, ...]:
        """Fetch & cache SeqRecords via Entrez.efetch.

        Args:
            ids: Tuple of record IDs to fetch.
            db: Entrez database name.
            rettype: Return format (e.g., 'gb', 'fasta').
            parser: SeqIO parser format.

        Returns:
            Tuple of SeqRecord objects.
        """
        if not ids:
            return ()
        with Entrez.efetch(
            db=db,
            id=",".join(ids),
            rettype=rettype,
            retmode="text",
        ) as handle:
            return tuple(SeqIO.parse(handle, parser))

    def _fetch_fasta_records(
        self, ids: Sequence[str], db: str
    ) -> List[SeqRecord]:
        """Wrapper to fetch FASTA via cached _efetch_records.

        Args:
            ids: Sequence of record IDs.
            db: Entrez database name.

        Returns:
            List of SeqRecord objects in FASTA format.
        """
        return list(self._efetch_records(tuple(ids), db, "fasta", "fasta"))

    def _fetch_genbank_records(
        self, ids: Sequence[str], db: str
    ) -> List[SeqRecord]:
        """Wrapper to fetch GenBank via cached _efetch_records.

        Args:
            ids: Sequence of record IDs.
            db: Entrez database name.

        Returns:
            List of SeqRecord objects in GenBank format.
        """
        return list(self._efetch_records(tuple(ids), db, "gb", "genbank"))

    @staticmethod
    @lru_cache(maxsize=None)
    def _fetch_linked_record_ids(
        gene_id: str, db: str, linkname: str
    ) -> List[str]:
        """Fetch linked record IDs via Entrez.elink().

        Args:
            gene_id: NCBI Gene ID.
            db: Target Entrez database.
            linkname: Link name for Entrez.elink.

        Returns:
            List of linked record IDs.
        """
        with Entrez.elink(
            dbfrom="gene", db=db, id=gene_id, linkname=linkname
        ) as handle:
            linksets = Entrez.read(handle)

        ids: List[str] = []
        for ls in linksets:
            for dbset in ls.get("LinkSetDb", []):
                for link in dbset.get("Link", []):
                    ids.append(link["Id"])
        return ids

    @staticmethod
    def _compute_exon_arrangement(
        exon_info: Dict[str, List[Dict[str, Any]]]
    ) -> Dict[str, Dict[str, List[int]]]:
        """Compute canonical and alternative exon‐number lists.

        Args:
            exon_info: Mapping from transcript ID to exon info list.

        Returns:
            Dict with 'canonical' and 'alternatives' exon arrangements.
        """
        if not exon_info:
            return {"canonical": {}, "alternatives": {}}

        # Compute total exon lengths per transcript
        lengths = {
            tid: sum(e["end"] - e["start"] + 1 for e in exons)
            for tid, exons in exon_info.items()
        }
        # The transcript with the longest total length is canonical
        canon = max(lengths, key=lengths.get)

        return {
            "canonical": {canon: [e["exon_number"] for e in exon_info[canon]]},
            "alternatives": {
                tid: [e["exon_number"] for e in exons]
                for tid, exons in exon_info.items()
                if tid != canon
            },
        }


# --- Example usage ---
"""
if __name__ == "__main__":

    # Initialize the fetcher (replace email with your own address if you like)
    fetcher = NCBIGeneSequenceFetcher(
        email="user@example.com",
        organism="Homo sapiens",
        assembly="GRCh38.p14",
    )

    gene_name = "CYP2C19"

    try:
        # Fetch raw SeqRecord objects and exon info
        raw_data = fetcher.fetch_all_ncbi_sequences(gene_name)

        # Genomic
        print("=== Genomic FASTA record (first 100 nt) ===")
        print(raw_data["genomic"].seq[:100], "\n")

        # CCDS (NCBI doesn’t provide CCDS via this fetcher)
        print("=== CCDS records ===")
        if raw_data.get("ccds"):
            for ccds_id, ccds_seq in raw_data["ccds"].items():
                print(f"{ccds_id}: {ccds_seq[:50]}... (length {len(ccds_seq)})")
        else:
            print("No CCDS records found.\n")

        # cDNA transcripts
        print("=== Transcript cDNA records ===")
        for rec in raw_data["cdna"]:
            print(f"{rec.id}: {str(rec.seq)[:50]}... (length {len(rec.seq)})")
        print()

        # Protein translations
        print("=== Transcript Protein records ===")
        for rec in raw_data["protein"]:
            if rec.seq:
                print(f"{rec.id}: {str(rec.seq)[:50]}... (length {len(rec.seq)})")
            else:
                print(f"{rec.id}: No protein sequence available.")
        print()

        # MANE Select (if available)
        if raw_data["mane_cdna"]:
            print("=== MANE Select cDNA ===")
            print(f"{raw_data['mane_cdna'].id}: {str(raw_data['mane_cdna'].seq)[:50]}... "
                  f"(length {len(raw_data['mane_cdna'].seq)})\n")
        else:
            print("No MANE Select cDNA record found.\n")

        if raw_data["mane_protein"]:
            print("=== MANE Select protein ===")
            print(f"{raw_data['mane_protein'].id}: {str(raw_data['mane_protein'].seq)[:50]}... "
                  f"(length {len(raw_data['mane_protein'].seq)})\n")
        else:
            print("No MANE Select protein record found.\n")

        # Exon information
        print("=== Exon coordinates & sequences ===")
        for tid, exons in raw_data["exon_info"].items():
            print(f"{tid}: {len(exons)} exons")
            for exon in exons:
                print(f"  Exon {exon['exon_number']}: {exon['start']}-{exon['end']} "
                      f"(len {len(exon['sequence'])})")
            print()

        # Exon arrangements
        print("=== Exon arrangement ===")
        print("Canonical:", raw_data["exon_arrangements"]["canonical"])
        print("Alternatives:", raw_data["exon_arrangements"]["alternatives"])

    except Exception as e:
        print(f"Error fetching data for {gene_name}: {e}")

    # — or — get everything as plain strings & dicts:
    print("\n--- Plain Python dict via fetch_all() ---")
    all_plain = fetcher.fetch_all(gene_name)
    print("Genomic seq length:", len(all_plain["genomic_sequence"]))
    print("cDNA transcripts:", list(all_plain["cdna"].keys()))
    print("Protein IDs:", list(all_plain["protein"].keys()))
    print("MANE cDNA:", all_plain["mane_cdna"])
    print("MANE protein:", all_plain["mane_protein"])
    print("Exon info keys:", all_plain["exon_info"].keys())
"""
