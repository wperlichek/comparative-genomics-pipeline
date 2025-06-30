import pytest
from unittest.mock import patch
import logging

from comparative_genomics_pipeline.util.file_util import validate_genes_config


class TestConfigValidation:
    """Test suite for gene configuration validation functions."""
    
    def test_validate_genes_config_valid_minimal(self):
        """Test validation with minimal valid configuration."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "P35498"
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is True
    
    def test_validate_genes_config_valid_complete(self):
        """Test validation with complete valid configuration."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "P35498",
                    "entrez_protein_id": "NP_001165963.1"
                },
                {
                    "species": "Mus musculus",
                    "uniprot_id": "Q99MZ9",
                    "entrez_protein_id": "NP_001074730.1"
                }
            ],
            "DEPDC5": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "O75140"
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is True
    
    def test_validate_genes_config_not_dict(self):
        """Test validation with non-dictionary input."""
        invalid_configs = [
            "not a dict",
            ["not", "a", "dict"],
            123,
            None
        ]
        
        for config in invalid_configs:
            result = validate_genes_config(config)
            assert result is False
    
    def test_validate_genes_config_empty_dict(self):
        """Test validation with empty configuration."""
        config = {}
        
        result = validate_genes_config(config)
        assert result is False
    
    def test_validate_genes_config_invalid_gene_names(self):
        """Test validation with invalid gene names."""
        invalid_configs = [
            {
                "": [{"species": "Homo sapiens", "uniprot_id": "P35498"}]  # Empty gene name
            },
            {
                "   ": [{"species": "Homo sapiens", "uniprot_id": "P35498"}]  # Whitespace gene name
            },
            {
                123: [{"species": "Homo sapiens", "uniprot_id": "P35498"}]  # Non-string gene name
            }
        ]
        
        for config in invalid_configs:
            result = validate_genes_config(config)
            assert result is False
    
    def test_validate_genes_config_ortholog_list_not_array(self):
        """Test validation when ortholog list is not an array."""
        invalid_configs = [
            {
                "SCN1A": "not an array"
            },
            {
                "SCN1A": {"species": "Homo sapiens"}  # Object instead of array
            },
            {
                "SCN1A": 123  # Number instead of array
            }
        ]
        
        for config in invalid_configs:
            result = validate_genes_config(config)
            assert result is False
    
    def test_validate_genes_config_empty_ortholog_list(self):
        """Test validation with empty ortholog list."""
        config = {
            "SCN1A": []
        }
        
        result = validate_genes_config(config)
        assert result is False
    
    def test_validate_genes_config_ortholog_not_object(self):
        """Test validation when ortholog entries are not objects."""
        invalid_configs = [
            {
                "SCN1A": ["not an object"]
            },
            {
                "SCN1A": [123]
            },
            {
                "SCN1A": [None]
            },
            {
                "SCN1A": [
                    {"species": "Homo sapiens", "uniprot_id": "P35498"},
                    "invalid ortholog"  # Mixed valid and invalid
                ]
            }
        ]
        
        for config in invalid_configs:
            result = validate_genes_config(config)
            assert result is False
    
    def test_validate_genes_config_missing_species(self):
        """Test validation with missing species field."""
        config = {
            "SCN1A": [
                {
                    "uniprot_id": "P35498"
                    # Missing species field
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is False
    
    def test_validate_genes_config_invalid_species(self):
        """Test validation with invalid species values."""
        invalid_configs = [
            {
                "SCN1A": [{"species": "", "uniprot_id": "P35498"}]  # Empty species
            },
            {
                "SCN1A": [{"species": "   ", "uniprot_id": "P35498"}]  # Whitespace species
            },
            {
                "SCN1A": [{"species": 123, "uniprot_id": "P35498"}]  # Non-string species
            },
            {
                "SCN1A": [{"species": None, "uniprot_id": "P35498"}]  # None species
            }
        ]
        
        for config in invalid_configs:
            result = validate_genes_config(config)
            assert result is False
    
    def test_validate_genes_config_missing_both_ids(self):
        """Test validation when both uniprot_id and entrez_protein_id are missing."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens"
                    # Missing both uniprot_id and entrez_protein_id
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is False
    
    def test_validate_genes_config_empty_ids(self):
        """Test validation when IDs are empty strings."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "",
                    "entrez_protein_id": ""
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is False
    
    def test_validate_genes_config_whitespace_ids(self):
        """Test validation when IDs are whitespace."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "   ",
                    "entrez_protein_id": "   "
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is False
    
    def test_validate_genes_config_uniprot_id_only_valid(self):
        """Test validation with only uniprot_id (should be valid)."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "P35498"
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is True
    
    def test_validate_genes_config_entrez_id_only_valid(self):
        """Test validation with only entrez_protein_id (should be valid)."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens",
                    "entrez_protein_id": "NP_001165963.1"
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is True
    
    def test_validate_genes_config_mixed_valid_invalid(self):
        """Test validation with mix of valid and invalid entries."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "P35498"
                },
                {
                    "species": "",  # Invalid: empty species
                    "uniprot_id": "Q99MZ9"
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is False
    
    def test_validate_genes_config_multiple_genes_valid(self):
        """Test validation with multiple valid genes."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "P35498"
                },
                {
                    "species": "Mus musculus",
                    "uniprot_id": "Q99MZ9"
                }
            ],
            "DEPDC5": [
                {
                    "species": "Homo sapiens",
                    "entrez_protein_id": "NP_055883.2"
                }
            ],
            "KCNQ2": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "O43526",
                    "entrez_protein_id": "NP_742105.1"
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is True
    
    def test_validate_genes_config_one_gene_invalid_fails_all(self):
        """Test that one invalid gene fails entire validation."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "P35498"
                }
            ],
            "INVALID_GENE": [
                {
                    "species": ""  # Invalid: empty species
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is False
    
    def test_validate_genes_config_realistic_epilepsy_genes(self):
        """Test validation with realistic epilepsy gene configuration."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "P35498",
                    "entrez_protein_id": "NP_001165963.1"
                },
                {
                    "species": "Pan troglodytes",
                    "uniprot_id": "H2QM44"
                },
                {
                    "species": "Mus musculus",
                    "uniprot_id": "Q99MZ9",
                    "entrez_protein_id": "NP_001074730.1"
                },
                {
                    "species": "Rattus norvegicus",
                    "uniprot_id": "P63088"
                },
                {
                    "species": "Danio rerio",
                    "entrez_protein_id": "NP_571641.1"
                }
            ],
            "DEPDC5": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "O75140",
                    "entrez_protein_id": "NP_055883.2"
                },
                {
                    "species": "Pan troglodytes",
                    "uniprot_id": "H2QNC8"
                },
                {
                    "species": "Mus musculus",
                    "uniprot_id": "Q6P9R2"
                },
                {
                    "species": "Rattus norvegicus",
                    "uniprot_id": "Q5U2N0"
                },
                {
                    "species": "Danio rerio",
                    "uniprot_id": "A0A0R4IQF6"
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is True
    
    @patch('comparative_genomics_pipeline.util.file_util.logger')
    def test_validate_genes_config_logging_success(self, mock_logger):
        """Test that successful validation logs appropriate message."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "P35498"
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is True
        
        # Check that success message was logged
        mock_logger.info.assert_called_with("Configuration validation passed: 1 genes configured")
    
    @patch('comparative_genomics_pipeline.util.file_util.logger')
    def test_validate_genes_config_logging_failure(self, mock_logger):
        """Test that failed validation logs error messages."""
        config = {
            "": [  # Invalid gene name
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "P35498"
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is False
        
        # Check that error messages were logged
        mock_logger.error.assert_called()
        error_calls = [call[0][0] for call in mock_logger.error.call_args_list]
        assert "Configuration validation failed" in error_calls
    
    def test_validate_genes_config_extra_fields_allowed(self):
        """Test that extra fields in ortholog objects are allowed."""
        config = {
            "SCN1A": [
                {
                    "species": "Homo sapiens",
                    "uniprot_id": "P35498",
                    "extra_field": "extra_value",
                    "another_field": 123
                }
            ]
        }
        
        result = validate_genes_config(config)
        assert result is True  # Extra fields should be allowed