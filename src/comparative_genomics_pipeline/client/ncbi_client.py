import httpx
import logging

logger = logging.getLogger(__name__)

class NCBIClient:
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    def __init__(self):
        self.client = httpx.AsyncClient()

    async def fetch_protein_fasta_by_entrez_id(self, entrez_protein_id: str = "") -> str:
        if not entrez_protein_id:
            logger.warning("Must provide Entrez Protein ID to fetch FASTA sequence")
            return ""

        params = {
            "db": "protein",
            "id": entrez_protein_id,
            "rettype": "fasta",
            "retmode": "text",
        }

        try:
            response = await self.client.get(self.BASE_URL, params=params)
            response.raise_for_status()
            return response.text
        except httpx.HTTPStatusError as e:
            logger.error(f"HTTP error for {response.url}: {e.response.status_code} {e.response.text}")
            return ""
        except httpx.RequestError as e:
            logger.error(f"Network error while requesting {self.BASE_URL}: {e}")
            return ""

    async def close(self):
        await self.client.aclose()