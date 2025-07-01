import httpx
import logging
import csv
from pathlib import Path
from typing import Optional

from .s3_client import S3Client
from ..config.aws_config import get_aws_config

logger = logging.getLogger(__name__)


class UniProtClient:
    BASE_URL = "https://rest.uniprot.org/uniprotkb/"

    def __init__(self):
        self.client = httpx.AsyncClient()
        
        # Initialize S3 caching (optional)
        try:
            aws_config = get_aws_config()
            self.s3_client = S3Client(aws_config)
            logger.info("S3 caching enabled for UniProt sequences")
        except ValueError as e:
            logger.info(f"S3 caching disabled: {e}")
            self.s3_client = None

    async def fetch_protein_fasta_sequence_by_accession_id(
        self, accession_id: str = ""
    ) -> str:

        if accession_id == "":
            logger.warning(
                "Must provide protein's accession id to get its FASTA sequence"
            )
            return ""
        
        # Check S3 cache first
        if self.s3_client:
            cached_sequence = await self.s3_client.get_sequence(accession_id)
            if cached_sequence:
                logger.debug(f"Using cached sequence for {accession_id}")
                return cached_sequence
        
        # Fallback to UniProt API
        url = f"{self.BASE_URL}{accession_id}.fasta"
        try:
            response = await self.client.get(url)
            response.raise_for_status()
            fasta_sequence = response.text
            
            # Cache the result if S3 is available
            if self.s3_client and fasta_sequence:
                await self.s3_client.put_sequence(accession_id, fasta_sequence)
                
            return fasta_sequence
            
        except httpx.HTTPStatusError as e:
            logger.error(
                f"HTTP error for {url}: {e.response.status_code} {e.response.text}"
            )
            return ""
        except httpx.RequestError as e:
            logger.error(f"Network error while requesting {url}: {e}")
            return ""

    async def fetch_protein_variants_by_accession_id(
        self, accession_id: str, output_dir: str, gene_name: str = None
    ):
        url = f"{self.BASE_URL}{accession_id}.json"
        logger.info(f"Requesting UniProt JSON for {accession_id}: {url}")
        try:
            response = await self.client.get(url)
            logger.info(f"HTTP status for {accession_id}: {response.status_code}")
            response.raise_for_status()
            data = response.json()
            logger.info(f"Top-level keys: {list(data.keys())}")
        except Exception as e:
            logger.error(f"Failed to fetch variants for {accession_id}: {e}")
            return

        features = data.get("features", [])
        logger.info(f"Total features: {len(features)}")
        types = set(feat.get("type", "") for feat in features)
        logger.info(f"Unique feature types: {types}")

        variants = []
        for i, feat in enumerate(features):
            logger.debug(f"Feature {i}: {feat}")
            if "variant" in feat.get("type", "").lower():
                pos = feat.get("location", {}).get("start", "")
                orig = feat.get("wildType", "")
                var = feat.get("alternativeSequence", "")
                desc = feat.get("description", "")
                variants.append(
                    {
                        "position": pos,
                        "original": orig,
                        "variant": var,
                        "description": desc,
                    }
                )
        logger.info(f"Extracted {len(variants)} variant features for {accession_id}")
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        filename = f"{gene_name}_{accession_id}_variants.csv" if gene_name else f"{accession_id}_variants.csv"
        out_path = Path(output_dir) / filename
        with open(out_path, "w", newline="") as csvfile:
            writer = csv.DictWriter(
                csvfile, fieldnames=["position", "original", "variant", "description"]
            )
            writer.writeheader()
            for v in variants:
                writer.writerow(v)
        identifier = gene_name if gene_name else accession_id
        logger.info(f"Saved {len(variants)} variants for {identifier} to {out_path}")

    async def close(self):
        await self.client.aclose()
