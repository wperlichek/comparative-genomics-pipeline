import pytest
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock, call
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade

from comparative_genomics_pipeline.service.biopython_service import visualize_and_save_trees


class TestPhylogeneticAnalysis:
    """Test suite for phylogenetic tree processing functions."""
    
    def create_test_newick_file(self, newick_string, filename="test_tree.nwk"):
        """Helper method to create test Newick tree files."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.nwk') as tmp:
            tmp.write(newick_string)
            return Path(tmp.name)
    
    def create_test_tree_object(self):
        """Helper method to create a test Bio.Phylo Tree object."""
        # Create a simple tree: ((A,B)C,D)E
        root = Clade(name="E")
        internal = Clade(name="C")
        internal.clades.extend([Clade(name="A"), Clade(name="B")])
        root.clades.extend([internal, Clade(name="D")])
        return Tree(root=root)
    
    def test_visualize_and_save_trees_single_valid_tree(self):
        """Test successful visualization of a single valid tree."""
        newick = "((A:0.1,B:0.2):0.05,C:0.3);"
        tree_file = self.create_test_newick_file(newick)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            
            try:
                result = visualize_and_save_trees([tree_file], output_dir)
                assert result is True
                
                # Check that PNG file was created
                png_file = output_dir / f"{tree_file.stem}.png"
                assert png_file.exists()
                
            finally:
                tree_file.unlink()
    
    def test_visualize_and_save_trees_multiple_trees(self):
        """Test batch processing of multiple tree files."""
        newick1 = "((A:0.1,B:0.2):0.05,C:0.3);"
        newick2 = "((D:0.1,E:0.2):0.05,F:0.3);"
        
        tree1 = self.create_test_newick_file(newick1, "tree1.nwk")
        tree2 = self.create_test_newick_file(newick2, "tree2.nwk")
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            
            try:
                result = visualize_and_save_trees([tree1, tree2], output_dir)
                assert result is True
                
                # Check that PNG files were created (names based on temp file stems)
                png_files = list(output_dir.glob("*.png"))
                assert len(png_files) == 2  # Two files should be created
                
            finally:
                tree1.unlink()
                tree2.unlink()
    
    def test_visualize_and_save_trees_auto_discover_files(self):
        """Test automatic discovery of tree files in directory."""
        newick = "((A:0.1,B:0.2):0.05,C:0.3);"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            
            # Create tree file directly in output directory
            tree_file = output_dir / "auto_tree.nwk"
            tree_file.write_text(newick)
            
            # Call without specifying tree_files (should auto-discover)
            result = visualize_and_save_trees(tree_files=None, output_dir=output_dir)
            assert result is True
            
            # Check that PNG file was created
            png_file = output_dir / "auto_tree.png"
            assert png_file.exists()
    
    def test_visualize_and_save_trees_invalid_newick(self):
        """Test error handling for malformed Newick files."""
        # Use truly invalid Newick that BioPython will reject
        invalid_newick = "((("  # Incomplete parentheses
        tree_file = self.create_test_newick_file(invalid_newick)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            
            try:
                result = visualize_and_save_trees([tree_file], output_dir)
                # Should return False because no trees were successfully processed
                assert result is False
                
            finally:
                tree_file.unlink()
    
    def test_visualize_and_save_trees_nonexistent_file(self):
        """Test error handling for non-existent tree files."""
        nonexistent_file = Path("/tmp/nonexistent_tree.nwk")
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            
            result = visualize_and_save_trees([nonexistent_file], output_dir)
            # Should return False because no trees were successfully processed
            assert result is False
    
    def test_visualize_and_save_trees_empty_tree_list(self):
        """Test behavior with empty tree file list."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            
            result = visualize_and_save_trees([], output_dir)
            assert result is True  # Should succeed with no files to process
    
    def test_visualize_and_save_trees_no_files_found(self):
        """Test behavior when no .nwk files are found in directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            
            # Create a non-.nwk file
            (output_dir / "not_a_tree.txt").write_text("some text")
            
            result = visualize_and_save_trees(tree_files=None, output_dir=output_dir)
            assert result is True  # Should succeed with warning about no files
    
    def test_visualize_and_save_trees_permission_error(self):
        """Test error handling for output directory permission issues."""
        newick = "((A:0.1,B:0.2):0.05,C:0.3);"
        tree_file = self.create_test_newick_file(newick)
        
        # Use a directory that doesn't exist and can't be created
        output_dir = Path("/root/inaccessible_dir")
        
        try:
            result = visualize_and_save_trees([tree_file], output_dir)
            assert result is False  # Should fail due to permission error
            
        finally:
            tree_file.unlink()
    
    @patch('comparative_genomics_pipeline.service.biopython_service.Phylo.read')
    @patch('comparative_genomics_pipeline.service.biopython_service.Phylo.draw')
    @patch('comparative_genomics_pipeline.service.biopython_service.plt.savefig')
    def test_visualize_and_save_trees_matplotlib_error(self, mock_savefig, mock_draw, mock_read):
        """Test error handling for matplotlib plotting errors."""
        newick = "((A:0.1,B:0.2):0.05,C:0.3);"
        tree_file = self.create_test_newick_file(newick)
        
        # Mock successful tree reading but failed plotting
        mock_read.return_value = self.create_test_tree_object()
        mock_draw.side_effect = Exception("Matplotlib error")
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            
            try:
                result = visualize_and_save_trees([tree_file], output_dir)
                # Should return False because plotting failed
                assert result is False
                
                # Should have attempted to read the tree
                mock_read.assert_called_once()
                
            finally:
                tree_file.unlink()
    
    @patch('comparative_genomics_pipeline.service.biopython_service.plt.close')
    def test_visualize_and_save_trees_figure_cleanup(self, mock_close):
        """Test that matplotlib figures are properly closed to prevent memory leaks."""
        newick = "((A:0.1,B:0.2):0.05,C:0.3);"
        tree_file = self.create_test_newick_file(newick)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            
            try:
                result = visualize_and_save_trees([tree_file], output_dir)
                assert result is True
                
                # Should have called plt.close() to cleanup figure
                mock_close.assert_called()
                
            finally:
                tree_file.unlink()
    
    def test_visualize_and_save_trees_complex_newick(self):
        """Test with complex Newick tree including branch lengths and bootstrap values."""
        # Complex tree with branch lengths and bootstrap support
        complex_newick = "((A:0.1,B:0.2)0.95:0.05,(C:0.15,D:0.25)0.80:0.1)0.99:0.0;"
        tree_file = self.create_test_newick_file(complex_newick)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            
            try:
                result = visualize_and_save_trees([tree_file], output_dir)
                assert result is True
                
                png_file = output_dir / f"{tree_file.stem}.png"
                assert png_file.exists()
                assert png_file.stat().st_size > 0  # File should not be empty
                
            finally:
                tree_file.unlink()
    
    def test_visualize_and_save_trees_species_names(self):
        """Test with realistic species names in tree."""
        # Tree with species names similar to those in the pipeline
        species_newick = "((Homo_sapiens:0.1,Pan_troglodytes:0.2):0.05,Mus_musculus:0.3);"
        tree_file = self.create_test_newick_file(species_newick)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            
            try:
                result = visualize_and_save_trees([tree_file], output_dir)
                assert result is True
                
                png_file = output_dir / f"{tree_file.stem}.png"
                assert png_file.exists()
                
            finally:
                tree_file.unlink()
    
    @patch('comparative_genomics_pipeline.service.biopython_service.path_config')
    def test_visualize_and_save_trees_default_paths(self, mock_path_config):
        """Test default path configuration usage."""
        with tempfile.TemporaryDirectory() as tmpdir:
            mock_path_config.TREES_OUTPUT_DIR = Path(tmpdir)
            
            newick = "((A:0.1,B:0.2):0.05,C:0.3);"
            tree_file = Path(tmpdir) / "default_tree.nwk"
            tree_file.write_text(newick)
            
            # Call with default parameters
            result = visualize_and_save_trees()
            assert result is True
            
            png_file = Path(tmpdir) / "default_tree.png"
            assert png_file.exists()
    
    def test_visualize_and_save_trees_edge_case_single_node(self):
        """Test tree with single node (edge case)."""
        single_node_newick = "A;"
        tree_file = self.create_test_newick_file(single_node_newick)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            
            try:
                result = visualize_and_save_trees([tree_file], output_dir)
                assert result is True
                
                png_file = output_dir / f"{tree_file.stem}.png"
                assert png_file.exists()
                
            finally:
                tree_file.unlink()