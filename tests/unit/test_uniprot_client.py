import pytest
from unittest.mock import AsyncMock, patch, MagicMock
import httpx
from comparative_genomics_pipeline.client.uniprot_client import UniProtClient


class TestUniProtClient:
    
    @pytest.fixture
    def client(self):
        """Create UniProt client instance with S3 disabled."""
        with patch('comparative_genomics_pipeline.client.uniprot_client.get_aws_config', side_effect=ValueError("Test mode - no AWS")):
            return UniProtClient()
    
    @pytest.fixture
    def mock_response(self):
        """Mock httpx response."""
        response = AsyncMock()
        response.status_code = 200
        response.text = ">sp|P35498|SCN1A_HUMAN Test sequence\nMATEST"
        response.json.return_value = {
            "features": [
                {
                    "type": "Natural variant",
                    "location": {"start": 100},
                    "wildType": "A",
                    "alternativeSequence": "T",
                    "description": "Test variant"
                }
            ]
        }
        return response

    @pytest.mark.unit
    async def test_fetch_protein_fasta_sequence_success(self, client):
        """Test successful FASTA sequence retrieval."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = ">sp|P35498|SCN1A_HUMAN Test sequence\nMATEST"
        mock_response.raise_for_status = MagicMock()
        
        with patch.object(client.client, 'get', return_value=mock_response):
            result = await client.fetch_protein_fasta_sequence_by_accession_id("P35498")
            assert result == ">sp|P35498|SCN1A_HUMAN Test sequence\nMATEST"
            client.client.get.assert_called_once_with(
                "https://rest.uniprot.org/uniprotkb/P35498.fasta"
            )

    @pytest.mark.unit
    async def test_fetch_protein_fasta_sequence_empty_id(self, client):
        """Test handling of empty accession ID."""
        result = await client.fetch_protein_fasta_sequence_by_accession_id("")
        assert result == ""

    @pytest.mark.unit
    async def test_fetch_protein_fasta_sequence_http_error(self, client):
        """Test handling of HTTP errors."""
        # Use MagicMock for synchronous methods
        mock_response = MagicMock()
        mock_response.status_code = 404
        mock_response.text = "Not found"
        
        # Create a proper HTTPStatusError with the mock response
        http_error = httpx.HTTPStatusError(
            message="Not found", request=MagicMock(), response=mock_response
        )
        mock_response.raise_for_status.side_effect = http_error
        
        with patch.object(client.client, 'get', return_value=mock_response):
            result = await client.fetch_protein_fasta_sequence_by_accession_id("INVALID")
            assert result == ""

    @pytest.mark.unit
    async def test_fetch_protein_fasta_sequence_network_error(self, client):
        """Test handling of network errors."""
        with patch.object(client.client, 'get', side_effect=httpx.RequestError("Network error")):
            result = await client.fetch_protein_fasta_sequence_by_accession_id("P35498")
            assert result == ""

    @pytest.mark.unit
    async def test_fetch_protein_variants_success(self, client, temp_data_dir):
        """Test successful variant retrieval."""
        # Create a properly configured mock response for variants
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.raise_for_status = MagicMock()  # No exception raised
        mock_response.json.return_value = {
            "features": [
                {
                    "type": "Natural variant",
                    "location": {"start": 100},
                    "wildType": "A",
                    "alternativeSequence": "T",
                    "description": "Test variant"
                }
            ]
        }
        
        with patch.object(client.client, 'get', return_value=mock_response):
            await client.fetch_protein_variants_by_accession_id("P35498", str(temp_data_dir))
            
            # Check that CSV file was created
            csv_file = temp_data_dir / "P35498_variants.csv"
            assert csv_file.exists()
            
            # Check CSV content
            content = csv_file.read_text()
            assert "position,original,variant,description" in content
            assert "100,A,T,Test variant" in content

    @pytest.mark.unit
    async def test_fetch_protein_variants_network_error(self, client, temp_data_dir):
        """Test handling of network errors in variant fetching."""
        with patch.object(client.client, 'get', side_effect=httpx.RequestError("Network error")):
            # Should not raise exception, just log error
            await client.fetch_protein_variants_by_accession_id("P35498", str(temp_data_dir))

    @pytest.mark.unit
    async def test_close(self, client):
        """Test client cleanup."""
        with patch.object(client.client, 'aclose') as mock_close:
            await client.close()
            mock_close.assert_called_once()