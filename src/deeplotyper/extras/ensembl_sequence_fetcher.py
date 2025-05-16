import logging
from typing import Any, Dict, List, Optional, Tuple

import requests
from requests import Session
from requests.exceptions import HTTPError

logger = logging.getLogger(__name__)


class EnsemblGeneDataFetcher:
    """Fetch gene and transcript data from the Ensembl REST API."""

    _BASE_URL: str = "https://rest.ensembl.org"
    _HEADERS: Dict[str, str] = {"Content-Type": "application/json"}

    def __init__(self, gene_symbol: str, species: str = "human") -> None:
        """Initialize fetcher with gene symbol and species.

        Args:
            gene_symbol: Gene symbol (e.g., "CYP2C19").
            species: Species name (e.g., "human" or "mus_musculus").
        """
        self.gene_symbol = gene_symbol
        self.species = species
        self._session: Session = requests.Session()
        self._session.headers.update(self._HEADERS)

    def _get_json(
        self, endpoint: str, params: Optional[Dict[str, Any]] = None
    ) -> Any:
        """GET a JSON payload from an Ensembl endpoint.

        Args:
            endpoint: Relative path (e.g., "/lookup/symbol/...").
            params: Query parameters.

        Returns:
            Parsed JSON response.

        Raises:
            HTTPError: If the request failed.
        """
        url = f"{self._BASE_URL}{endpoint}"
        response = self._session.get(url, params=params)
        response.raise_for_status()
        return response.json()

    def _fetch_sequence(
        self, path: str, params: Optional[Dict[str, Any]] = None
    ) -> str:
        """Fetch a sequence (cDNA, protein, or genomic) as a string.

        Args:
            path: Endpoint path for the sequence.
            params: Query parameters (e.g., {"type": "cdna"}).

        Returns:
            Sequence string, or an empty string if not found.
        """
        data = self._get_json(path, params=params)
        return data.get("seq", "")

    @staticmethod
    def _normalize_genomic_info (raw: Dict[str, Any]) -> Dict[str, Any]:
        """Map Ensembl lookup JSON into our unified genomic_info schema."""
        return {
            "assembly": raw.get("assembly_name"),
            "annotation_release": None,
            "seq_region_name": raw.get("seq_region_name"),
            "seq_region_accession": None,
            "start": raw.get("start"),
            "end": raw.get("end")
        }

    def get_genomic_sequence(self) -> Tuple[Dict[str, Any], str]:
        """Lookup gene coordinates and fetch its genomic region as FASTA.

        Returns:
            A tuple of:
                gene_lookup_json: JSON dict from the lookup endpoint.
                genomic_sequence: FASTA sequence string.
        """
        lookup = self._get_json(
            f"/lookup/symbol/{self.species}/{self.gene_symbol}"
        )
        region = (
            f"{lookup['seq_region_name']}:"
            f"{lookup['start']}..{lookup['end']}:"
            f"{lookup['strand']}"
        )
        sequence = self._fetch_sequence(
            f"/sequence/region/{self.species}/{region}"
        )
        return lookup, sequence

    def get_ccds_records(self) -> Dict[str, str]:
        """Retrieve CCDS records via Ensembl xrefs.

        Returns:
            Mapping of CCDS ID to cleaned cDNA sequence
            (trailing 'N's stripped). Skips non-CCDS or failed fetches.
        """
        try:
            xrefs = self._get_json(
                f"/xrefs/symbol/{self.species}/{self.gene_symbol}",
                params={"external_db": "CCDS"},
            )
        except HTTPError:
            return {}

        records: Dict[str, str] = {}
        for entry in xrefs:
            ccds_id = entry.get("id", "")
            if not ccds_id.startswith("CCDS"):
                logger.debug("Skipping non-CCDS xref: %r", ccds_id)
                continue

            try:
                seq = self._fetch_sequence(
                    f"/sequence/id/{ccds_id}", params={"type": "cdna"}
                )
            except HTTPError as err:
                logger.warning(
                    "Failed to fetch CCDS sequence for %s: %s", ccds_id, err
                )
                continue

            records[ccds_id] = seq.rstrip("N")
        return records

    def get_transcript_info(
        self
    ) -> Dict[str, Dict[str, Dict[str, Optional[str]]]]:
        """Fetch transcripts, their cDNA & protein, and pick canonical one.

        Returns:
            Dict with keys 'canonical' and 'alternatives', each mapping
            transcript ID to a dict with keys 'cdna', 'protein', 'biotype'.

        Raises:
            ValueError: If no transcripts are found for the gene.
        """
        gene_data = self._get_json(
            f"/lookup/symbol/{self.species}/{self.gene_symbol}",
            params={"expand": 1},
        )
        transcripts = gene_data.get("Transcript", [])
        if not transcripts:
            raise ValueError(
                f"No transcripts found for gene '{self.gene_symbol}'."
            )

        info: Dict[str, Dict[str, Optional[str]]] = {}
        for tx in transcripts:
            tx_id = tx["id"]
            biotype = tx.get("biotype", "unknown")
            cdna_seq = self._fetch_sequence(
                f"/sequence/id/{tx_id}", params={"type": "cdna"}
            )

            protein_seq: Optional[str] = None
            translation = tx.get("Translation")
            if translation and translation.get("id"):
                protein_seq = self._fetch_sequence(
                    f"/sequence/id/{translation['id']}",
                    params={"type": "protein"},
                )

            info[tx_id] = {
                "cdna": cdna_seq,
                "protein": protein_seq,
                "biotype": biotype,
            }

        canonical_id = next(
            (t["id"] for t in transcripts if t.get("is_canonical") == 1),
            max(info, key=lambda tid: len(info[tid]["cdna"])),
        )
        canonical = {canonical_id: info.pop(canonical_id)}
        return {"canonical": canonical, "alternatives": info}

    @staticmethod
    def combine_transcript_info(
        transcript_info: Dict[str, Dict[str, Dict[str, Optional[str]]]]
    ) -> Tuple[Dict[str, str], Dict[str, Optional[str]]]:
        """Flatten cDNA and protein sequences for all transcripts.

        Args:
            transcript_info: Dict returned by get_transcript_info().

        Returns:
            A tuple of:
                cdna_map: Mapping from transcript ID to cDNA string.
                protein_map: Mapping from transcript ID to protein string or None.
        """
        cdna_map: Dict[str, str] = {}
        protein_map: Dict[str, Optional[str]] = {}
        for group in ("canonical", "alternatives"):
            for tx_id, attrs in transcript_info[group].items():
                cdna_map[tx_id] = attrs["cdna"]
                protein_map[tx_id] = attrs["protein"]
        return cdna_map, protein_map

    def get_exon_info(
            self, transcript_ids: Optional[List[str]] = None
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Fetch exon features and sequences for transcripts, assigning per-transcript exon ranks.

        Args:
            transcript_ids: List of transcript IDs. If None, uses the canonical one.

        Returns:
            Mapping of transcript ID to a list of exon dicts with keys:
                exon_id, exon_number (rank), seq_region, start, end, strand, sequence.
        """
        if transcript_ids is None:
            transcript_ids = list(self.get_transcript_info()["canonical"].keys())

        exon_info: Dict[str, List[Dict[str, Any]]] = {}
        for tx_id in transcript_ids:
            # Fetch the transcript with its nested Exon children:
            tx_data = self._get_json(f"/lookup/id/{tx_id}", params={"expand": 1})
            exons = tx_data.get("Exon", [])

            # Sort by genomic start (so rank = position in transcript)
            sorted_exons = sorted(exons, key=lambda e: e["start"])

            entries: List[Dict[str, Any]] = []
            for exon_number, exon in enumerate(sorted_exons, start=1):
                region = (
                    f"{exon['seq_region_name']}:"
                    f"{exon['start']}..{exon['end']}:"
                    f"{exon['strand']}"
                )
                seq = self._fetch_sequence(f"/sequence/region/{self.species}/{region}")

                entries.append({
                    "exon_id": exon["id"],
                    "exon_number": exon_number,
                    "seq_region": exon["seq_region_name"],
                    "start": exon["start"],
                    "end": exon["end"],
                    "strand": exon["strand"],
                    "sequence": seq,
                })

            exon_info[tx_id] = entries

        return exon_info

    def get_exon_arrangement_ensembl(
        self
    ) -> Dict[str, Dict[str, List[int]]]:
        """Return true exon-rank lists for canonical vs. alternative transcripts.

        Uses each transcript’s own exon_number (rank) rather than positional index.
        """
        tx_info = self.get_transcript_info()
        all_ids = list(tx_info["canonical"]) + list(tx_info["alternatives"])
        exon_data = self.get_exon_info(all_ids)

        return {
            "canonical": {
                tx_id: [exon["exon_number"] for exon in exon_data[tx_id]]
                for tx_id in tx_info["canonical"]
            },
            "alternatives": {
                tx_id: [exon["exon_number"] for exon in exon_data[tx_id]]
                for tx_id in tx_info["alternatives"]
            },
        }

    def fetch_all_ensembl_sequences(self) -> Dict[str, Any]:
        """Retrieve comprehensive sequence data for the gene.

        Returns:
            Dict containing:
                gene_record: FASTA genomic sequence string,
                genomic_info: lookup JSON dict,
                ccds: CCDS ID→sequence map,
                cdna: transcript ID→cDNA map,
                protein: transcript ID→protein map,
                exon_arrangements: exon-number mappings,
                exon_info_canonical: exon info for the canonical transcript.
        """
        gene_lookup, genome_seq = self.get_genomic_sequence()
        ccds = self.get_ccds_records()
        tx_info = self.get_transcript_info()
        cdna_map, protein_map = self.combine_transcript_info(tx_info)
        arrangements = self.get_exon_arrangement_ensembl()
        exon_info_canonical = self.get_exon_info(
            list(tx_info["canonical"])
        )

        return {
            "gene_record": genome_seq,
            "cdna": cdna_map,
            "ccds": ccds,
            "protein": protein_map,
            "exon_arrangements": arrangements,
            "exon_info_canonical": exon_info_canonical,
            "genomic_info": gene_lookup,
        }

    def fetch_all(self) -> Dict[str, Any]:
        # get the raw dict
        raw = self.fetch_all_ensembl_sequences()

        # Ensembl returns genome seq as string, lookup JSON in raw['genomic_info']
        # cdna/protein are already maps {tx:seq}, and ccds is a map too
        return {
            "genomic_sequence": raw["gene_record"],
            "genomic_info":      self._normalize_genomic_info(raw["genomic_info"]),
            "cdna":               raw["cdna"],
            "protein":            raw["protein"],
            "ccds":               raw["ccds"],
            # Ensembl doesn’t have MANE, so we fill with None
            "mane_cdna":          None,
            "mane_protein":       None,
            # fallback: Ensembl only returns exon_info for canonical
            # if you want *all* transcripts you could call get_exon_info on all ids
            "exon_info":          raw["exon_info_canonical"],
            "exon_arrangements":  raw["exon_arrangements"],
        }

# --- Example usage ---
"""
if __name__ == "__main__":

    gene_symbol = "CYP2C19"
    fetcher = EnsemblGeneDataFetcher(gene_symbol)

    try:
        gene_data = fetcher.get_gene_data()
        print("=== Genomic FASTA record (first 100 nt) ===")
        print(gene_data["genomic"][:100], "\n")

        print("=== CCDS records ===")
        if gene_data["ccds"]:
            for ccds_id, ccds_seq in gene_data["ccds"].items():
                print(f"{ccds_id}: {ccds_seq[:50]}... (length {len(ccds_seq)})")
        else:
            print("No CCDS records found.")

        print("\n=== Transcript cDNA records ===")
        for tid, seq_cdna in gene_data["cdna"].items():
            print(f"{tid}: {seq_cdna[:50]}... (length {len(seq_cdna)})")

        print("\n=== Transcript Protein records ===")
        for tid, seq_protein in gene_data["protein"].items():
            if seq_protein:
                print(f"{tid}: {seq_protein[:50]}... (length {len(seq_protein)})")
            else:
                print(f"{tid}: No protein sequence available.")
    except Exception as e:
        print(f"Error: {e}")


    cdna_records = gene_data["cdna"]
    protein_translations = fetcher.translate_cdna_records(cdna_records)
    print(protein_translations)
"""

