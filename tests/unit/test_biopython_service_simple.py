import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
import tempfile
from comparative_genomics_pipeline.service.biopython_service import (
    visualize_and_save_trees
)


class TestBiopythonServiceSimple:
    """Simplified tests for biopython service that actually work."""

    @pytest.mark.unit
    @patch('comparative_genomics_pipeline.service.biopython_service.Phylo.read')
    @patch('comparative_genomics_pipeline.service.biopython_service.Phylo.draw')
    @patch('matplotlib.pyplot.figure')
    @patch('matplotlib.pyplot.savefig')
    @patch('matplotlib.pyplot.close')
    @patch('builtins.print')
    def test_visualize_and_save_trees_with_files(self, mock_print, mock_close, mock_savefig, mock_figure, mock_draw, mock_phylo_read):
        """Test tree visualization with provided tree files."""
        # Setup mocks
        mock_tree = MagicMock()
        mock_phylo_read.return_value = mock_tree
        mock_fig = MagicMock()
        mock_figure.return_value = mock_fig
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            tree_file = temp_path / "test_tree.nwk"
            tree_file.touch()  # Create empty file
            
            visualize_and_save_trees([tree_file], temp_path)
            
            # Verify tree was read and processed
            mock_phylo_read.assert_called_once_with(tree_file, "newick")
            mock_figure.assert_called_once_with(figsize=(10, 5))
            mock_draw.assert_called_once_with(mock_tree, do_show=False)
            mock_savefig.assert_called_once_with(temp_path / "test_tree.png", bbox_inches="tight")
            mock_close.assert_called_once_with(mock_fig)

    @pytest.mark.unit 
    def test_module_imports(self):
        """Test that the module imports correctly."""
        from comparative_genomics_pipeline.service import biopython_service
        assert hasattr(biopython_service, 'visualize_and_save_trees')
        assert hasattr(biopython_service, 'compute_conservation_scores')