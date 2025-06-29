import pytest
import tempfile
import shutil
from pathlib import Path
from typing import Dict, Any
from unittest.mock import AsyncMock

@pytest.fixture
def temp_data_dir():
    """Create a temporary directory for test data."""
    temp_dir = tempfile.mkdtemp()
    yield Path(temp_dir)
    shutil.rmtree(temp_dir)

@pytest.fixture
def sample_gene_config() -> Dict[str, Any]:
    """Sample gene configuration for testing."""
    return {
        "SCN1A": [
            {
                "species": "Homo sapiens",
                "uniprot_id": "P35498",
                "entrez_protein_id": "ENSP05155037736.1"
            },
            {
                "species": "Mus musculus", 
                "uniprot_id": "A2APX8",
                "entrez_protein_id": "ENSMUSP00215039336.1"
            }
        ]
    }

@pytest.fixture
def mock_httpx_client():
    """Mock httpx async client."""
    return AsyncMock()

@pytest.fixture
def sample_fasta_content():
    """Sample FASTA content for testing."""
    return """>sp|P35498|SCN1A_HUMAN Sodium channel protein type 1 subunit alpha
MAASDSEYRTRSEAETLSITDMEAGTDVQKADGDFVQGQHQEVSKVQGTGTDSGAFQHGPQATP
>sp|A2APX8|SCN1A_MOUSE Sodium channel protein type 1 subunit alpha  
MAASDSEYRTRSEAETLSITDMEAGTDVQKADGDFVQGQHQEVSKVQGTGTDSGAFQHGPQATP"""

@pytest.fixture 
def sample_alignment_data():
    """Sample alignment data for conservation analysis."""
    return {
        "sequences": ["MAASD", "MAASD", "MAAAD"],
        "positions": list(range(5)),
        "species": ["Human", "Mouse", "Chicken"]
    }