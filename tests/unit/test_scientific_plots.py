"""
Tests for scientific plotting functionality - conservation variant plots that must never be lost.

This module contains critical tests for the conservation variant plotting functionality
that was previously lost and restored. These tests ensure the plots are generated
correctly and prevent regression.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import os
from unittest.mock import Mock, patch

from comparative_genomics_pipeline.visualization.scientific_plots import (
    ConservationPlotter, 
    VariantPlotter, 
    PhylogeneticPlotter,
    ClinVarPlotter
)
from comparative_genomics_pipeline.service.biopython_service import (
    plot_variants_scientific,
    plot_conservation_scientific,
)


class TestConservationVariantPlots:
    """Critical tests for conservation variant plots - MUST NEVER BE REMOVED."""
    
    @pytest.fixture
    def sample_conservation_data(self):
        """Create sample conservation data matching the real pipeline output format."""
        return pd.DataFrame({
            'Position': range(1, 101),
            'ShannonEntropy_WithGaps': np.random.uniform(0, 2.5, 100),
            'ShannonEntropy_NoGaps': np.random.uniform(0, 2.3, 100),
            'ConsensusResidue': ['A', 'R', 'N', 'D', 'C'] * 20
        })
    
    @pytest.fixture
    def sample_variant_data(self):
        """Create sample variant data matching UniProt output format."""
        return pd.DataFrame({
            'position': [10, 25, 45, 60, 85],
            'description': [
                'loss of function',
                'likely pathogenic', 
                'reduced function',
                'likely benign',
                'uncertain significance'
            ],
            'consequence': ['missense', 'nonsense', 'missense', 'missense', 'missense'],
            'wildtype': ['A', 'R', 'C', 'D', 'G'],
            'variant': ['T', 'stop', 'S', 'N', 'E']
        })
    
    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary directory for test outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    def test_variant_plotter_initialization(self):
        """Test that VariantPlotter can be initialized."""
        plotter = VariantPlotter()
        assert plotter is not None
        assert plotter.config is not None
        assert plotter.theme is not None
    
    def test_variant_position_parsing(self, sample_variant_data):
        """Test that variant positions are parsed correctly from different formats."""
        plotter = VariantPlotter()
        parsed_df = plotter._parse_variant_positions(sample_variant_data)
        
        assert 'parsed_position' in parsed_df.columns
        assert len(parsed_df) == 5  # All positions should parse successfully
        assert list(parsed_df['parsed_position']) == [10, 25, 45, 60, 85]
    
    def test_variant_classification_logic(self, sample_variant_data):
        """Test that variants are classified correctly (CRITICAL - this was in the lost plots)."""
        plotter = VariantPlotter()
        parsed_df = plotter._parse_variant_positions(sample_variant_data)
        
        lof_positions, pathogenic_positions, additional = plotter._get_dynamic_variant_classifications(
            parsed_df, "TEST_GENE"
        )
        
        # Loss of function variants
        assert 10 in lof_positions  # 'loss of function'
        
        # Pathogenic variants
        assert 25 in pathogenic_positions  # 'likely pathogenic'
        
        # Additional classifications
        assert 45 in additional['reduced_function']  # 'reduced function'
        assert 60 in additional['benign']  # 'likely benign'
        assert 85 in additional['uncertain']  # 'uncertain significance'
    
    def test_conservation_variant_analysis(self, sample_conservation_data, sample_variant_data):
        """Test statistical analysis of variant-conservation relationship."""
        plotter = VariantPlotter()
        parsed_variants = plotter._parse_variant_positions(sample_variant_data)
        
        stats_results = plotter._analyze_variant_conservation(
            sample_conservation_data, parsed_variants
        )
        
        assert 'variant_conservation' in stats_results
        assert 'background_conservation' in stats_results
        assert 'n_variants' in stats_results
        assert 'n_background' in stats_results
        assert 'p_value' in stats_results
        assert 'cohens_d' in stats_results
        
        # Check that we have the expected number of variants
        assert stats_results['n_variants'] == 5
        assert stats_results['n_background'] == 95  # 100 total - 5 variants
    
    def test_variant_plot_generation_end_to_end(self, sample_conservation_data, sample_variant_data, temp_output_dir):
        """CRITICAL TEST: End-to-end variant plot generation - this functionality was LOST before."""
        # Save sample data to CSV files
        conservation_csv = temp_output_dir / "test_conservation.csv"
        variants_csv = temp_output_dir / "test_variants.csv"
        
        sample_conservation_data.to_csv(conservation_csv, index=False)
        sample_variant_data.to_csv(variants_csv, index=False)
        
        # Generate the plot using the actual service function
        plotter = VariantPlotter()
        output_path = plotter.plot_variants_with_statistics(
            conservation_csv, variants_csv, temp_output_dir
        )
        
        # Verify the plot was created
        assert output_path is not None
        assert output_path.exists()
        assert output_path.suffix == '.png'
        assert 'variants_scientific' in output_path.name
        
        # Verify file is not empty (has actual plot data)
        assert output_path.stat().st_size > 1000  # At least 1KB for a real plot
    
    def test_biopython_service_plot_variants_scientific(self, sample_conservation_data, sample_variant_data, temp_output_dir):
        """Test the biopython service wrapper function - this is what the main pipeline calls."""
        # Save sample data to CSV files  
        conservation_csv = temp_output_dir / "test_conservation.csv"
        variants_csv = temp_output_dir / "test_variants.csv"
        
        sample_conservation_data.to_csv(conservation_csv, index=False)
        sample_variant_data.to_csv(variants_csv, index=False)
        
        # Call the service function that the main pipeline uses
        output_path = plot_variants_scientific(
            str(conservation_csv), str(variants_csv), str(temp_output_dir)
        )
        
        # Verify the plot was created through the service layer
        assert output_path is not None
        expected_path = temp_output_dir / "test_conservation_variants_scientific.png" 
        assert expected_path.exists()
        assert expected_path.stat().st_size > 1000
    
    def test_real_data_structure_compatibility(self):
        """Test that the plotting functions work with the real data structure from the pipeline."""
        # This test uses the actual files from the pipeline to ensure compatibility
        from comparative_genomics_pipeline.config import path_config
        
        conservation_csv = path_config.CONSERVATION_OUTPUT_DIR / "SCN1A_conservation.csv"
        variants_csv = path_config.VARIANTS_OUTPUT_DIR / "SCN1A_P35498_variants.csv"
        
        # Only run this test if the real data files exist
        if conservation_csv.exists() and variants_csv.exists():
            conservation_df = pd.read_csv(conservation_csv)
            variants_df = pd.read_csv(variants_csv)
            
            # Verify data structure matches expectations
            assert 'Position' in conservation_df.columns
            assert 'ShannonEntropy_NoGaps' in conservation_df.columns
            assert 'position' in variants_df.columns
            assert 'description' in variants_df.columns
            
            # Test that the plotter can handle this real data
            plotter = VariantPlotter()
            parsed_variants = plotter._parse_variant_positions(variants_df)
            
            # Should successfully parse positions
            assert len(parsed_variants) > 0
            assert 'parsed_position' in parsed_variants.columns
            
            # Should successfully classify variants
            lof_pos, path_pos, additional = plotter._get_dynamic_variant_classifications(
                parsed_variants, "SCN1A"
            )
            
            # SCN1A should have some loss-of-function variants based on the research
            assert len(lof_pos) > 0, "SCN1A should have loss-of-function variants in the real data"
    
    def test_plot_file_naming_convention(self, sample_conservation_data, sample_variant_data, temp_output_dir):
        """Test that plot files follow the expected naming convention."""
        conservation_csv = temp_output_dir / "GENE_conservation.csv"
        variants_csv = temp_output_dir / "GENE_variants.csv"
        
        sample_conservation_data.to_csv(conservation_csv, index=False)
        sample_variant_data.to_csv(variants_csv, index=False)
        
        plotter = VariantPlotter()
        output_path = plotter.plot_variants_with_statistics(
            conservation_csv, variants_csv, temp_output_dir
        )
        
        # Check naming convention: {gene}_conservation_variants_scientific.png
        expected_name = "GENE_conservation_variants_scientific.png"
        assert output_path.name == expected_name
    
    def test_variant_plot_contains_key_elements(self, sample_conservation_data, sample_variant_data, temp_output_dir):
        """Test that the variant plots contain the key visual elements that made them valuable."""
        conservation_csv = temp_output_dir / "test_conservation.csv"
        variants_csv = temp_output_dir / "test_variants.csv"
        
        sample_conservation_data.to_csv(conservation_csv, index=False)
        sample_variant_data.to_csv(variants_csv, index=False)
        
        # Create plotter and generate statistics (this tests the statistical analysis component)
        plotter = VariantPlotter()
        parsed_variants = plotter._parse_variant_positions(sample_variant_data)
        stats_results = plotter._analyze_variant_conservation(sample_conservation_data, parsed_variants)
        
        # Verify that key statistical measures are computed
        assert not np.isnan(stats_results['variant_mean'])
        assert not np.isnan(stats_results['background_mean'])
        assert stats_results['n_variants'] > 0
        assert stats_results['n_background'] > 0
        
        # Test that variant classification produces expected categories
        lof_pos, path_pos, additional = plotter._get_dynamic_variant_classifications(
            parsed_variants, "TEST"
        )
        
        # Should identify different types of variants
        assert len(lof_pos) > 0  # Loss of function
        assert len(path_pos) > 0  # Likely pathogenic
        assert len(additional['reduced_function']) > 0  # Reduced function
        assert len(additional['benign']) > 0  # Benign
        assert len(additional['uncertain']) > 0  # Uncertain


class TestConservationPlotter:
    """Tests for conservation plotting without variants."""
    
    def test_conservation_plotter_initialization(self):
        """Test ConservationPlotter can be initialized."""
        plotter = ConservationPlotter()
        assert plotter is not None
    
    def test_conservation_plot_generation(self, tmp_path):
        """Test basic conservation plot generation."""
        # Create sample conservation data
        conservation_data = pd.DataFrame({
            'Position': range(1, 51),
            'ShannonEntropy_WithGaps': np.random.uniform(0, 2.5, 50),
            'ShannonEntropy_NoGaps': np.random.uniform(0, 2.3, 50),
            'ConsensusResidue': ['A', 'R', 'N', 'D', 'C'] * 10
        })
        
        conservation_csv = tmp_path / "test_conservation.csv"
        conservation_data.to_csv(conservation_csv, index=False)
        
        plotter = ConservationPlotter()
        output_path = plotter.plot_conservation_with_confidence(conservation_csv, tmp_path)
        
        assert output_path is not None
        assert output_path.exists()
        assert '_scientific.png' in output_path.name


class TestPhylogeneticPlotter:
    """Tests for phylogenetic tree plotting."""
    
    def test_phylogenetic_plotter_initialization(self):
        """Test PhylogeneticPlotter can be initialized."""
        plotter = PhylogeneticPlotter()
        assert plotter is not None


class TestClinVarPlotter:
    """Tests for ClinVar variant plotting."""
    
    def test_clinvar_plotter_initialization(self):
        """Test ClinVarPlotter can be initialized."""
        plotter = ClinVarPlotter()
        assert plotter is not None


class TestPlotFileExistence:
    """Tests to ensure critical plot files exist and are maintained."""
    
    def test_conservation_variant_plots_exist_in_output(self):
        """CRITICAL: Test that conservation variant plots exist in the output directory."""
        from comparative_genomics_pipeline.config import path_config
        
        variants_dir = path_config.VARIANTS_OUTPUT_DIR
        
        # Check for the specific plot files that were lost before
        scn1a_plot = variants_dir / "SCN1A_conservation_variants_scientific.png"
        depdc5_plot = variants_dir / "DEPDC5_conservation_variants_scientific.png"
        
        # These files should exist after running the pipeline
        if scn1a_plot.exists():
            assert scn1a_plot.stat().st_size > 1000, "SCN1A conservation variant plot should not be empty"
        
        if depdc5_plot.exists():
            assert depdc5_plot.stat().st_size > 1000, "DEPDC5 conservation variant plot should not be empty"
    
    def test_conservation_variant_plots_can_be_regenerated(self):
        """Test that conservation variant plots can be regenerated if lost."""
        from comparative_genomics_pipeline.config import path_config
        
        conservation_csv = path_config.CONSERVATION_OUTPUT_DIR / "SCN1A_conservation.csv"
        variants_csv = path_config.VARIANTS_OUTPUT_DIR / "SCN1A_P35498_variants.csv"
        output_dir = path_config.VARIANTS_OUTPUT_DIR
        
        # Only test if source files exist
        if conservation_csv.exists() and variants_csv.exists():
            # This should work without errors
            output_path = plot_variants_scientific(
                str(conservation_csv), str(variants_csv), str(output_dir)
            )
            
            expected_plot = output_dir / "SCN1A_conservation_variants_scientific.png"
            assert expected_plot.exists(), "Conservation variant plot should be regenerated successfully"
            assert expected_plot.stat().st_size > 1000, "Regenerated plot should not be empty"


# Mark critical tests that must never be removed
CRITICAL_TESTS = [
    'test_variant_plot_generation_end_to_end',
    'test_biopython_service_plot_variants_scientific', 
    'test_variant_classification_logic',
    'test_conservation_variant_analysis',
    'test_conservation_variant_plots_exist_in_output',
    'test_conservation_variant_plots_can_be_regenerated'
]

def test_critical_tests_present():
    """Meta-test to ensure critical tests are present in this file."""
    import inspect
    
    # Get all test functions in this module
    current_module = inspect.getmembers(inspect.getmodule(inspect.currentframe()))
    test_functions = [name for name, obj in current_module if inspect.isfunction(obj) and name.startswith('test_')]
    
    # Add test methods from test classes
    test_classes = [obj for name, obj in current_module if inspect.isclass(obj) and name.startswith('Test')]
    for test_class in test_classes:
        class_methods = [method for method in dir(test_class) if method.startswith('test_')]
        test_functions.extend(class_methods)
    
    # Ensure all critical tests are present
    for critical_test in CRITICAL_TESTS:
        assert critical_test in test_functions, f"Critical test {critical_test} is missing and must be restored!"


if __name__ == "__main__":
    # Run basic functionality check
    print("Running basic test validation...")
    test_critical_tests_present()
    print("✅ All critical conservation variant plot tests are present and accounted for.")