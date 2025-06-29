import pytest
from unittest.mock import AsyncMock, patch, MagicMock
import httpx
from comparative_genomics_pipeline.client.ebi_client import EBIClient


class TestEBIClient:
    
    @pytest.fixture
    def client(self):
        """Create EBI client instance."""
        return EBIClient()
    
    @pytest.fixture
    def sample_fasta(self):
        """Sample FASTA sequence for testing."""
        return ">seq1\nMATEST\n>seq2\nMATEXT"

    @pytest.mark.unit
    async def test_submit_job_success(self, client, sample_fasta):
        """Test successful job submission."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "clustalo-test-job-123"
        mock_response.raise_for_status = MagicMock()
        
        with patch.object(client.client, 'post', return_value=mock_response):
            job_id = await client.submit_job(sample_fasta)
            assert job_id == "clustalo-test-job-123"
            client.client.post.assert_called_once_with(
                "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run/",
                data={"sequence": sample_fasta, "email": "author@gmail.com"}
            )

    @pytest.mark.unit
    async def test_submit_job_http_error(self, client, sample_fasta):
        """Test handling of HTTP errors during job submission."""
        mock_response = MagicMock()
        mock_response.status_code = 400
        mock_response.text = "Bad request"
        
        http_error = httpx.HTTPStatusError(
            message="Bad request", request=MagicMock(), response=mock_response
        )
        mock_response.raise_for_status.side_effect = http_error
        
        with patch.object(client.client, 'post', return_value=mock_response):
            job_id = await client.submit_job(sample_fasta)
            assert job_id is None

    @pytest.mark.unit
    async def test_submit_job_network_error(self, client, sample_fasta):
        """Test handling of network errors during job submission."""
        with patch.object(client.client, 'post', side_effect=httpx.RequestError("Network error")):
            job_id = await client.submit_job(sample_fasta)
            assert job_id is None

    @pytest.mark.unit
    async def test_check_status_success(self, client):
        """Test successful status check."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "FINISHED"
        mock_response.raise_for_status = MagicMock()
        
        with patch.object(client.client, 'get', return_value=mock_response):
            status = await client.check_status("test-job-123")
            assert status == "FINISHED"
            client.client.get.assert_called_once_with(
                "https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/test-job-123"
            )

    @pytest.mark.unit
    async def test_check_status_http_error(self, client):
        """Test handling of HTTP errors during status check."""
        mock_response = MagicMock()
        mock_response.status_code = 404
        mock_response.text = "Not found"
        
        http_error = httpx.HTTPStatusError(
            message="Not found", request=MagicMock(), response=mock_response
        )
        mock_response.raise_for_status.side_effect = http_error
        
        with patch.object(client.client, 'get', return_value=mock_response):
            status = await client.check_status("invalid-job")
            assert status is None

    @pytest.mark.unit
    async def test_get_result_success(self, client):
        """Test successful result retrieval."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = ">seq1\nMA-TEST\n>seq2\nMATEXT-"
        mock_response.raise_for_status = MagicMock()
        
        with patch.object(client.client, 'get', return_value=mock_response):
            result = await client.get_result("test-job-123", "fa")
            assert result == ">seq1\nMA-TEST\n>seq2\nMATEXT-"
            client.client.get.assert_called_once_with(
                "https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/test-job-123/fa"
            )

    @pytest.mark.unit
    async def test_get_result_http_error(self, client):
        """Test handling of HTTP errors during result retrieval."""
        mock_response = MagicMock()
        mock_response.status_code = 404
        mock_response.text = "Not found"
        
        http_error = httpx.HTTPStatusError(
            message="Not found", request=MagicMock(), response=mock_response
        )
        mock_response.raise_for_status.side_effect = http_error
        
        with patch.object(client.client, 'get', return_value=mock_response):
            result = await client.get_result("invalid-job", "fa")
            assert result is None

    @pytest.mark.unit
    async def test_get_phylogenetic_tree(self, client):
        """Test phylogenetic tree retrieval."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "(seq1:0.1,seq2:0.2);"
        mock_response.raise_for_status = MagicMock()
        
        with patch.object(client.client, 'get', return_value=mock_response):
            tree = await client.get_phylogenetic_tree("test-job-123")
            assert tree == "(seq1:0.1,seq2:0.2);"
            client.client.get.assert_called_once_with(
                "https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/test-job-123/phylotree"
            )

    @pytest.mark.unit
    async def test_close(self, client):
        """Test client cleanup."""
        with patch.object(client.client, 'aclose') as mock_close:
            await client.close()
            mock_close.assert_called_once()