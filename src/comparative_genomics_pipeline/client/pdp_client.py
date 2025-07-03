import httpx
from pathlib import Path
from ..config import path_config
import logging
import json
from typing import Dict, List, Any, Optional

logger = logging.getLogger(__name__)

class PDBClient:
    BASE_URL = "https://files.rcsb.org/download/"
    UNIPROT_API = "https://rest.uniprot.org/uniprotkb/"

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

    def fetch_uniprot_domains(self, accession: str) -> Optional[Dict[str, Any]]:
        """
        Fetch domain information from UniProt for a given accession.
        Returns domain data including positions and types.
        """
        url = f"{self.UNIPROT_API}{accession}.json"
        try:
            response = self.client.get(url)
            response.raise_for_status()
            data = response.json()
            
            domains = []
            features = data.get('features', [])
            
            for feature in features:
                # Include key structural features
                if feature.get('type') in ['Domain', 'Transmembrane', 'Repeat', 'Topological domain', 'Region', 'Intramembrane']:
                    location = feature.get('location', {})
                    start = location.get('start', {}).get('value')
                    end = location.get('end', {}).get('value')
                    
                    if start and end:
                        domains.append({
                            'type': feature.get('type'),
                            'description': feature.get('description', ''),
                            'start': start,
                            'end': end,
                            'length': end - start + 1
                        })
            
            protein_length = data.get('sequence', {}).get('length', 0)
            
            return {
                'accession': accession,
                'protein_length': protein_length,
                'domains': sorted(domains, key=lambda x: x['start'])
            }
            
        except httpx.HTTPStatusError as e:
            logger.error(f"HTTP error fetching UniProt data for {accession}: {e.response.status_code}")
        except httpx.RequestError as e:
            logger.error(f"Network error fetching UniProt data for {accession}: {e}")
        except Exception as e:
            logger.error(f"Error parsing UniProt data for {accession}: {e}")
        
        return None

    def save_domain_data(self, domain_data: Dict[str, Any], out_dir: Optional[Path] = None) -> Optional[Path]:
        """
        Save domain data to JSON file in structures directory.
        """
        if out_dir is None:
            out_dir = path_config.DATA_OUTPUT_DIR / "structures"
        out_dir.mkdir(parents=True, exist_ok=True)
        
        accession = domain_data.get('accession')
        if not accession:
            logger.error("No accession found in domain data")
            return None
            
        out_path = out_dir / f"{accession}_domains.json"
        
        try:
            with open(out_path, 'w') as f:
                json.dump(domain_data, f, indent=2)
            logger.info(f"Saved domain data for {accession} to {out_path}")
            return out_path
        except Exception as e:
            logger.error(f"Error saving domain data for {accession}: {e}")
            return None

    def close(self):
        self.client.close()
