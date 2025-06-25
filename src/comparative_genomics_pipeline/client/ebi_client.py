import httpx
import logging
from typing import Optional

logger = logging.getLogger(__name__)


class EBIClient:
    BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
    EMAIL = "williamperlichek@gmail.com"

    def __init__(self):
        self.client = httpx.AsyncClient()

    async def submit_job(self, fasta_sequence: str) -> Optional[str]:
        """
        Submit a FASTA sequence to Clustal Omega and return the job ID.
        """
        data = {"sequence": fasta_sequence, "email": self.EMAIL}
        try:
            response = await self.client.post(f"{self.BASE_URL}/run/", data=data)
            response.raise_for_status()
            job_id = response.text.strip()
            logger.info(f"Job submitted successfully: {job_id}")
            return job_id
        except httpx.HTTPStatusError as e:
            logger.error(
                f"HTTP error during job submission: {e.response.status_code} {e.response.text}"
            )
        except httpx.RequestError as e:
            logger.error(f"Network error during job submission: {e}")
        return None

    async def check_status(self, job_id: str) -> Optional[str]:
        """
        Check the status of a submitted Clustal Omega job.
        """
        try:
            response = await self.client.get(f"{self.BASE_URL}/status/{job_id}")
            response.raise_for_status()
            status = response.text.strip()
            logger.info(f"Job {job_id} status: {status}")
            return status
        except httpx.HTTPStatusError as e:
            logger.error(
                f"HTTP error checking status for {job_id}: {e.response.status_code} {e.response.text}"
            )
        except httpx.RequestError as e:
            logger.error(f"Network error checking status for {job_id}: {e}")
        return None

    async def get_result(self, job_id: str, result_type: str = "fa") -> Optional[str]:
        """
        Fetch a result from a completed Clustal Omega job.
        result_type: 'fa' for FASTA alignment, 'clustal' for CLUSTAL format, 'ph' for Newick tree, etc.
        """
        try:
            url = f"{self.BASE_URL}/result/{job_id}/{result_type}"
            response = await self.client.get(url)
            response.raise_for_status()
            return response.text
        except httpx.HTTPStatusError as e:
            logger.error(
                f"HTTP error fetching result {result_type} for {job_id}: {e.response.status_code} {e.response.text}"
            )
        except httpx.RequestError as e:
            logger.error(
                f"Network error fetching result {result_type} for {job_id}: {e}"
            )
        return None

    async def get_phylogenetic_tree(self, job_id: str) -> Optional[str]:
        """
        Fetch the phylogenetic tree in Newick format from the Clustal Omega job result.
        Returns the tree as a string, or None on error.
        """
        return await self.get_result(job_id, result_type="phylotree")

    async def close(self):
        await self.client.aclose()
