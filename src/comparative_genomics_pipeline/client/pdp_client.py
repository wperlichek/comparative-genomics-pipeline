import httpx
from pathlib import Path
from ..config import path_config
import logging

logger = logging.getLogger(__name__)

class PDBClient:
    BASE_URL = "https://files.rcsb.org/download/"

    def __init__(self):
        self.client = httpx.Client()

    def fetch_pdb(self, pdb_id: str, accession: str = None, out_dir: Path = None) -> Path:
        """
        Download a PDB file by its ID and save to the output/proteins directory.
        If accession is provided, include it in the filename (e.g., P35498_7DTD.pdb).
        Returns the path to the saved file, or None on failure.
        """
        if out_dir is None:
            out_dir = path_config.DATA_OUTPUT_DIR / "proteins"
        out_dir.mkdir(parents=True, exist_ok=True)
        url = f"{self.BASE_URL}{pdb_id}.pdb"
        if accession:
            out_path = out_dir / f"{accession}_{pdb_id}.pdb"
        else:
            out_path = out_dir / f"{pdb_id}.pdb"
        try:
            response = self.client.get(url)
            response.raise_for_status()
            with open(out_path, "w") as f:
                f.write(response.text)
            logger.info(f"Downloaded PDB {pdb_id} to {out_path}")
            return out_path
        except httpx.HTTPStatusError as e:
            logger.error(f"HTTP error fetching PDB {pdb_id}: {e.response.status_code} {e.response.text}")
        except httpx.RequestError as e:
            logger.error(f"Network error fetching PDB {pdb_id}: {e}")
        return None

    def close(self):
        self.client.close()
