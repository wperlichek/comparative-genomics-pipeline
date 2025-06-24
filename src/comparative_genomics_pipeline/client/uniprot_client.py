import httpx, logging

logger = logging.getLogger(__name__)


class UniProtClient:
    BASE_URL = "https://rest.uniprot.org/uniprotkb/"

    def __init__(self):
        self.client = httpx.AsyncClient()

    async def fetch_protein_fasta_sequence_by_accession_id(
        self, accession_id: str = ""
    ) -> str:

        if accession_id == "":
            logger.warning(
                "Must provide protein's accession id to get its FASTA sequence"
            )
            return ""
        else:
            url = f"{self.BASE_URL}{accession_id}"
            headers = {"Accept": "text/x-fasta"}
            try:
                response = await self.client.get(url, headers=headers)
                response.raise_for_status()
                return response.text
            except httpx.HTTPStatusError as e:
                logger.error(f"HTTP error for {url}: {e.response.status_code} {e.response.text}")
                return ""
            except httpx.RequestError as e:
                logger.error(f"Network error while requesting {url}: {e}")
                return ""

    async def close(self):
        await self.client.aclose()
