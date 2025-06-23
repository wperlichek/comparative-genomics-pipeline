import httpx

class UniProtClient:
    BASE_URL = "https://rest.uniprot.org/uniprotkb/"

    def __init__(self):
        self.client = httpx.AsyncClient()

    async def fetch_protein(self, query: str) -> str:
        url = f"{self.BASE_URL}search?query={query}&format=fasta"
        response = await self.client.get(url)
        response.raise_for_status()
        return response.text

    async def close(self):
        await self.client.aclose()