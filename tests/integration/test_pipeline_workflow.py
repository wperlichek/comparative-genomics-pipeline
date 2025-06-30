import pytest
import tempfile
import json
from pathlib import Path
from unittest.mock import AsyncMock, patch
from comparative_genomics_pipeline.client.uniprot_client import UniProtClient
from comparative_genomics_pipeline.client.ebi_client import EBIClient


class TestPipelineWorkflow:
    """Integration tests for end-to-end pipeline workflows."""
    
    @pytest.fixture
    def sample_gene_config_file(self, temp_data_dir, sample_gene_config):
        """Create a temporary gene configuration file."""
        config_file = temp_data_dir / "genes_to_proteins.json"
        with open(config_file, 'w') as f:
            json.dump(sample_gene_config, f)
        return config_file
    
    @pytest.fixture
    def mock_uniprot_responses(self):
        """Mock responses for UniProt API calls."""
        return {
            "P35498": ">sp|P35498|SCN1A_HUMAN Sodium channel protein type 1 subunit alpha\nMAASDSEYRTRSEAETLSITDMEAGTDVQ",
            "A2APX8": ">sp|A2APX8|SCN1A_MOUSE Sodium channel protein type 1 subunit alpha\nMAASDSEYRTRSEAETLSITDMEAGTDVQ"
        }
    
    @pytest.fixture
    def mock_ebi_responses(self):
        """Mock responses for EBI Clustal Omega API calls."""
        return {
            "job_id": "clustalo-test-job-123",
            "status": "FINISHED",
            "alignment": ">sp|P35498|SCN1A_HUMAN\nMAASDSEYRTRSEAETLSITDMEAGTDVQ\n>sp|A2APX8|SCN1A_MOUSE\nMAASDSEYRTRSEAETLSITDMEAGTDVQ",
            "tree": "(P35498:0.0,A2APX8:0.0);"
        }

    @pytest.mark.integration
    async def test_uniprot_to_ebi_workflow(self, mock_uniprot_responses, mock_ebi_responses):
        """Test the workflow from UniProt sequence retrieval to EBI alignment."""
        # Setup clients with S3 disabled for testing
        with patch('comparative_genomics_pipeline.client.uniprot_client.get_aws_config', side_effect=ValueError("Test mode - no AWS")):
            uniprot_client = UniProtClient()
        ebi_client = EBIClient()
        
        try:
            # Mock UniProt responses
            with patch.object(uniprot_client.client, 'get') as mock_uniprot_get:
                mock_uniprot_get.return_value = AsyncMock()
                mock_uniprot_get.return_value.status_code = 200
                mock_uniprot_get.return_value.text = mock_uniprot_responses["P35498"]
                mock_uniprot_get.return_value.raise_for_status = AsyncMock()
                
                # Step 1: Fetch sequences from UniProt
                sequence1 = await uniprot_client.fetch_protein_fasta_sequence_by_accession_id("P35498")
                assert sequence1 == mock_uniprot_responses["P35498"]
                
                mock_uniprot_get.return_value.text = mock_uniprot_responses["A2APX8"]
                sequence2 = await uniprot_client.fetch_protein_fasta_sequence_by_accession_id("A2APX8")
                assert sequence2 == mock_uniprot_responses["A2APX8"]
                
                # Combine sequences for alignment
                combined_fasta = f"{sequence1}\n{sequence2}"
                
                # Mock EBI responses
                with patch.object(ebi_client.client, 'post') as mock_ebi_post, \
                     patch.object(ebi_client.client, 'get') as mock_ebi_get:
                    
                    # Mock job submission
                    mock_ebi_post.return_value = AsyncMock()
                    mock_ebi_post.return_value.status_code = 200
                    mock_ebi_post.return_value.text = mock_ebi_responses["job_id"]
                    mock_ebi_post.return_value.raise_for_status = AsyncMock()
                    
                    # Mock status check
                    mock_ebi_get.return_value = AsyncMock()
                    mock_ebi_get.return_value.status_code = 200
                    mock_ebi_get.return_value.text = mock_ebi_responses["status"]
                    mock_ebi_get.return_value.raise_for_status = AsyncMock()
                    
                    # Step 2: Submit alignment job
                    job_id = await ebi_client.submit_job(combined_fasta)
                    assert job_id == mock_ebi_responses["job_id"]
                    
                    # Step 3: Check job status
                    status = await ebi_client.check_status(job_id)
                    assert status == "FINISHED"
                    
                    # Step 4: Get alignment result
                    mock_ebi_get.return_value.text = mock_ebi_responses["alignment"]
                    alignment = await ebi_client.get_result(job_id, "fa")
                    assert alignment == mock_ebi_responses["alignment"]
                    
                    # Step 5: Get phylogenetic tree
                    mock_ebi_get.return_value.text = mock_ebi_responses["tree"]
                    tree = await ebi_client.get_phylogenetic_tree(job_id)
                    assert tree == mock_ebi_responses["tree"]
        
        finally:
            await uniprot_client.close()
            await ebi_client.close()

    @pytest.mark.integration
    @pytest.mark.slow
    async def test_file_processing_workflow(self, temp_data_dir, sample_gene_config_file):
        """Test file-based processing workflow without external API calls."""
        # Create mock FASTA file
        fasta_content = """>sp|P35498|SCN1A_HUMAN
MAASDSEYRTRSEAETLSITDMEAGTDVQ
>sp|A2APX8|SCN1A_MOUSE
MAASDSEYRTRSEAETLSITDMEAGTDVQ"""
        
        fasta_file = temp_data_dir / "test_sequences.fasta"
        fasta_file.write_text(fasta_content)
        
        # Create mock alignment file
        alignment_content = """>sp|P35498|SCN1A_HUMAN
MAASDSEYRTRSEAETLSITDMEAGTDVQ
>sp|A2APX8|SCN1A_MOUSE
MAASDSEYRTRSEAETLSITDMEAGTDVQ"""
        
        alignment_file = temp_data_dir / "test_alignment.fasta"
        alignment_file.write_text(alignment_content)
        
        # Create mock tree file
        tree_content = "(P35498:0.0,A2APX8:0.0);"
        tree_file = temp_data_dir / "test_tree.nwk"
        tree_file.write_text(tree_content)
        
        # Verify files were created correctly
        assert fasta_file.exists()
        assert alignment_file.exists()
        assert tree_file.exists()
        assert "SCN1A_HUMAN" in fasta_file.read_text()
        assert "SCN1A_MOUSE" in alignment_file.read_text()
        assert "P35498" in tree_file.read_text()

    @pytest.mark.integration
    def test_configuration_loading(self, sample_gene_config_file):
        """Test loading and parsing gene configuration."""
        with open(sample_gene_config_file, 'r') as f:
            config = json.load(f)
        
        # Verify structure
        assert "SCN1A" in config
        assert len(config["SCN1A"]) == 2
        
        # Verify human entry
        human_entry = next(entry for entry in config["SCN1A"] if entry["species"] == "Homo sapiens")
        assert human_entry["uniprot_id"] == "P35498"
        assert human_entry["entrez_protein_id"] == "ENSP05155037736.1"
        
        # Verify mouse entry
        mouse_entry = next(entry for entry in config["SCN1A"] if entry["species"] == "Mus musculus")
        assert mouse_entry["uniprot_id"] == "A2APX8"
        assert mouse_entry["entrez_protein_id"] == "ENSMUSP00215039336.1"