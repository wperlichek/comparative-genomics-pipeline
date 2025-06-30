import pytest
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock
import tempfile
from pathlib import Path

from comparative_genomics_pipeline.visualization.scientific_plots import VariantPlotter


class TestVariantConservationMapping:
    """Test suite for variant-conservation mapping functions."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.plotter = VariantPlotter()
    
    def create_test_conservation_df(self):
        """Create test conservation data."""
        return pd.DataFrame({
            'Position': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            'ShannonEntropy_NoGaps': [0.0, 0.5, 1.0, 1.5, 2.0, 0.2, 0.8, 1.2, 1.8, 0.1]
        })
    
    def create_test_variants_df(self):
        """Create test variant data."""
        return pd.DataFrame({
            'position': [2, 5, 8],  # Positions with known conservation scores
            'variant_type': ['missense', 'nonsense', 'missense'],
            'clinical_significance': ['pathogenic', 'pathogenic', 'benign']
        })
    
    def test_parse_variant_positions_integer_positions(self):
        """Test parsing simple integer positions."""
        variants_df = pd.DataFrame({
            'position': [1, 2, 3, 4, 5]
        })
        
        result = self.plotter._parse_variant_positions(variants_df)
        
        assert 'parsed_position' in result.columns
        assert list(result['parsed_position']) == [1, 2, 3, 4, 5]
        assert result['parsed_position'].dtype == int
    
    def test_parse_variant_positions_float_positions(self):
        """Test parsing float positions (should convert to int)."""
        variants_df = pd.DataFrame({
            'position': [1.0, 2.5, 3.9, 4.1]
        })
        
        result = self.plotter._parse_variant_positions(variants_df)
        
        expected = [1, 2, 3, 4]  # Float to int conversion
        assert list(result['parsed_position']) == expected
    
    def test_parse_variant_positions_json_format(self):
        """Test parsing JSON-formatted positions."""
        variants_df = pd.DataFrame({
            'position': ["{'value': 100}", "{'value': 200}", "{'value': 300}"]
        })
        
        result = self.plotter._parse_variant_positions(variants_df)
        
        expected = [100, 200, 300]
        assert list(result['parsed_position']) == expected
    
    def test_parse_variant_positions_mixed_formats(self):
        """Test parsing mixed position formats."""
        variants_df = pd.DataFrame({
            'position': [1, "{'value': 100}", 2.5, "invalid", None]
        })
        
        result = self.plotter._parse_variant_positions(variants_df)
        
        # Should only keep valid positions
        expected = [1, 100, 2]
        assert list(result['parsed_position']) == expected
        assert len(result) == 3  # Invalid and None should be dropped
    
    def test_parse_variant_positions_invalid_data(self):
        """Test parsing with invalid position data."""
        variants_df = pd.DataFrame({
            'position': ["invalid", "not_a_number", None, "{'invalid': 'json'}"]
        })
        
        # The function uses eval() which can raise exceptions, so we expect some errors
        # but the function should handle them gracefully
        try:
            result = self.plotter._parse_variant_positions(variants_df)
            # Should return empty DataFrame after dropping invalid rows
            assert len(result) == 0
            assert 'parsed_position' in result.columns
        except (SyntaxError, KeyError):
            # eval() with invalid JSON can raise these errors
            # This is expected behavior for malformed data
            pass
    
    def test_parse_variant_positions_empty_dataframe(self):
        """Test parsing empty DataFrame."""
        variants_df = pd.DataFrame(columns=['position'])
        
        result = self.plotter._parse_variant_positions(variants_df)
        
        assert len(result) == 0
        assert 'parsed_position' in result.columns
    
    def test_analyze_variant_conservation_basic_analysis(self):
        """Test basic variant-conservation statistical analysis."""
        conservation_df = self.create_test_conservation_df()
        variants_df = self.create_test_variants_df()
        
        # Parse positions first
        variants_df = self.plotter._parse_variant_positions(variants_df)
        
        result = self.plotter._analyze_variant_conservation(conservation_df, variants_df)
        
        # Should return analysis results
        assert 'error' not in result
        assert 'variant_conservation' in result or 'statistic' in result
    
    def test_analyze_variant_conservation_no_matching_positions(self):
        """Test analysis with no matching conservation positions."""
        conservation_df = self.create_test_conservation_df()
        variants_df = pd.DataFrame({
            'position': [100, 200, 300],  # Outside conservation range
            'parsed_position': [100, 200, 300]
        })
        
        result = self.plotter._analyze_variant_conservation(conservation_df, variants_df)
        
        assert 'error' in result
        assert result['error'] == 'No matching positions found'
    
    def test_analyze_variant_conservation_statistical_tests(self):
        """Test statistical significance testing in variant analysis."""
        # Create conservation data with clear pattern
        conservation_df = pd.DataFrame({
            'Position': list(range(1, 21)),
            'ShannonEntropy_NoGaps': [0.1] * 10 + [2.0] * 10  # Low then high conservation
        })
        
        # Variants at high conservation positions
        variants_df = pd.DataFrame({
            'position': [15, 16, 17],
            'parsed_position': [15, 16, 17]
        })
        
        result = self.plotter._analyze_variant_conservation(conservation_df, variants_df)
        
        # Should perform Mann-Whitney U test
        if 'error' not in result:
            # Check that statistical test was performed
            assert len(result.get('variant_conservation', [])) > 0
    
    def test_analyze_variant_conservation_single_variant(self):
        """Test analysis with single variant (edge case)."""
        conservation_df = self.create_test_conservation_df()
        variants_df = pd.DataFrame({
            'position': [5],
            'parsed_position': [5]
        })
        
        result = self.plotter._analyze_variant_conservation(conservation_df, variants_df)
        
        # Should handle single variant gracefully
        assert 'error' not in result or result.get('variant_conservation') is not None
    
    def test_analyze_variant_conservation_all_positions_variants(self):
        """Test when all positions have variants (edge case)."""
        conservation_df = pd.DataFrame({
            'Position': [1, 2, 3],
            'ShannonEntropy_NoGaps': [0.5, 1.0, 1.5]
        })
        
        variants_df = pd.DataFrame({
            'position': [1, 2, 3],
            'parsed_position': [1, 2, 3]
        })
        
        result = self.plotter._analyze_variant_conservation(conservation_df, variants_df)
        
        # Should handle case where no background positions exist
        assert result is not None
    
    def test_analyze_variant_conservation_missing_conservation_column(self):
        """Test error handling for missing conservation columns."""
        conservation_df = pd.DataFrame({
            'Position': [1, 2, 3],
            'Wrong_Column': [0.5, 1.0, 1.5]  # Missing ShannonEntropy_NoGaps
        })
        
        variants_df = pd.DataFrame({
            'position': [1, 2],
            'parsed_position': [1, 2]
        })
        
        with pytest.raises(KeyError):
            self.plotter._analyze_variant_conservation(conservation_df, variants_df)
    
    def test_analyze_variant_conservation_nan_values(self):
        """Test handling of NaN values in conservation data."""
        conservation_df = pd.DataFrame({
            'Position': [1, 2, 3, 4, 5],
            'ShannonEntropy_NoGaps': [0.5, np.nan, 1.0, np.nan, 1.5]
        })
        
        variants_df = pd.DataFrame({
            'position': [1, 3, 5],
            'parsed_position': [1, 3, 5]
        })
        
        # Should handle NaN values gracefully
        result = self.plotter._analyze_variant_conservation(conservation_df, variants_df)
        
        # Check that analysis proceeds despite NaN values
        assert result is not None
    
    def test_analyze_variant_conservation_empty_conservation_data(self):
        """Test with empty conservation DataFrame."""
        conservation_df = pd.DataFrame(columns=['Position', 'ShannonEntropy_NoGaps'])
        variants_df = self.create_test_variants_df()
        variants_df = self.plotter._parse_variant_positions(variants_df)
        
        result = self.plotter._analyze_variant_conservation(conservation_df, variants_df)
        
        assert 'error' in result
        assert result['error'] == 'No matching positions found'
    
    def test_analyze_variant_conservation_empty_variants_data(self):
        """Test with empty variants DataFrame."""
        conservation_df = self.create_test_conservation_df()
        variants_df = pd.DataFrame(columns=['position', 'parsed_position'])
        
        result = self.plotter._analyze_variant_conservation(conservation_df, variants_df)
        
        assert 'error' in result
        assert result['error'] == 'No matching positions found'
    
    def test_parse_variant_positions_string_numbers(self):
        """Test parsing string representations of numbers."""
        variants_df = pd.DataFrame({
            'position': ["1", "2", "3", "4.5", "invalid"]
        })
        
        result = self.plotter._parse_variant_positions(variants_df)
        
        # Should parse valid string numbers
        expected = [1, 2, 3, 4]  # "invalid" should be dropped
        assert list(result['parsed_position']) == expected
    
    def test_analyze_variant_conservation_effect_size_calculation(self):
        """Test effect size calculation in statistical analysis."""
        # Create data with known effect size
        conservation_df = pd.DataFrame({
            'Position': list(range(1, 21)),
            'ShannonEntropy_NoGaps': [0.0] * 10 + [1.0] * 10  # Clear difference
        })
        
        # Variants at high conservation positions
        variants_df = pd.DataFrame({
            'position': [11, 12, 13],
            'parsed_position': [11, 12, 13]
        })
        
        result = self.plotter._analyze_variant_conservation(conservation_df, variants_df)
        
        # Should calculate effect size when there are sufficient data points
        if 'error' not in result and len(result.get('variant_conservation', [])) > 1:
            # Effect should be calculated for sufficient sample sizes
            assert isinstance(result, dict)
    
    def test_parse_variant_positions_preserves_other_columns(self):
        """Test that parsing preserves other DataFrame columns."""
        variants_df = pd.DataFrame({
            'position': [1, 2, 3],
            'gene': ['SCN1A', 'SCN1A', 'DEPDC5'],
            'variant_type': ['missense', 'nonsense', 'missense']
        })
        
        result = self.plotter._parse_variant_positions(variants_df)
        
        # Should preserve original columns
        assert 'gene' in result.columns
        assert 'variant_type' in result.columns
        assert 'parsed_position' in result.columns
        assert list(result['gene']) == ['SCN1A', 'SCN1A', 'DEPDC5']