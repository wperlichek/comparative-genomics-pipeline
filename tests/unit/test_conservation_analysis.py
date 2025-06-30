import pytest
import tempfile
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from comparative_genomics_pipeline.service.biopython_service import compute_conservation_scores, compute_conservation_for_all_msas


class TestConservationAnalysis:
    """Test suite for Shannon entropy conservation analysis functions."""
    
    def create_test_alignment_file(self, sequences, filename="test_alignment.fasta"):
        """Helper method to create test MSA files."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
            for i, seq in enumerate(sequences):
                tmp.write(f">seq_{i}\n{seq}\n")
            return Path(tmp.name)
    
    def test_compute_conservation_scores_perfect_conservation(self):
        """Test Shannon entropy calculation for perfectly conserved positions."""
        # Create alignment with perfectly conserved positions
        sequences = ["AAAA", "AAAA", "AAAA"]
        alignment_file = self.create_test_alignment_file(sequences)
        
        try:
            result_file = compute_conservation_scores(alignment_file)
            assert result_file is not None
            assert result_file.exists()
            
            # Read and validate results
            import pandas as pd
            df = pd.read_csv(result_file)
            
            # Perfect conservation should have entropy = 0
            assert len(df) == 4  # 4 positions
            assert all(df["ShannonEntropy_WithGaps"] == 0.0)
            assert all(df["ShannonEntropy_NoGaps"] == 0.0)
            assert list(df["Position"]) == [1, 2, 3, 4]
            
        finally:
            alignment_file.unlink()
            if result_file and result_file.exists():
                result_file.unlink()
    
    def test_compute_conservation_scores_maximum_diversity(self):
        """Test Shannon entropy for maximally diverse positions."""
        # Create alignment with 4 different amino acids per position
        sequences = ["AAAA", "CCCC", "GGGG", "TTTT"]
        alignment_file = self.create_test_alignment_file(sequences)
        
        try:
            result_file = compute_conservation_scores(alignment_file)
            assert result_file is not None
            
            import pandas as pd
            df = pd.read_csv(result_file)
            
            # Maximum entropy for 4 equal frequencies = log2(4) = 2.0
            expected_entropy = 2.0
            assert len(df) == 4
            assert all(abs(df["ShannonEntropy_WithGaps"] - expected_entropy) < 0.001)
            assert all(abs(df["ShannonEntropy_NoGaps"] - expected_entropy) < 0.001)
            
        finally:
            alignment_file.unlink()
            if result_file and result_file.exists():
                result_file.unlink()
    
    def test_compute_conservation_scores_with_gaps(self):
        """Test entropy calculation handling gaps correctly."""
        # Position 1: A,A,A,- (should differ with/without gaps)
        # Position 2: A,A,A,A (should be same with/without gaps)
        sequences = ["AA", "AA", "AA", "-A"]
        alignment_file = self.create_test_alignment_file(sequences)
        
        try:
            result_file = compute_conservation_scores(alignment_file)
            assert result_file is not None
            
            import pandas as pd
            df = pd.read_csv(result_file)
            
            assert len(df) == 2
            
            # Position 1: with gaps has A(3) and -(1), without gaps has A(3)
            pos1 = df[df["Position"] == 1].iloc[0]
            assert pos1["ShannonEntropy_WithGaps"] > 0  # Should have some entropy
            assert pos1["ShannonEntropy_NoGaps"] == 0   # Perfect conservation without gaps
            
            # Position 2: should be identical (all A's)
            pos2 = df[df["Position"] == 2].iloc[0]
            assert pos2["ShannonEntropy_WithGaps"] == 0
            assert pos2["ShannonEntropy_NoGaps"] == 0
            
        finally:
            alignment_file.unlink()
            if result_file and result_file.exists():
                result_file.unlink()
    
    def test_compute_conservation_scores_known_entropy_values(self):
        """Test against manually calculated Shannon entropy values."""
        # Position with known distribution: A(2), C(1), G(1) = 4 sequences
        # Expected entropy = -(2/4 * log2(2/4) + 1/4 * log2(1/4) + 1/4 * log2(1/4))
        # = -(0.5 * (-1) + 0.25 * (-2) + 0.25 * (-2)) = -(-0.5 - 0.5 - 0.5) = 1.5
        sequences = ["A", "A", "C", "G"]
        alignment_file = self.create_test_alignment_file(sequences)
        
        try:
            result_file = compute_conservation_scores(alignment_file)
            assert result_file is not None
            
            import pandas as pd
            df = pd.read_csv(result_file)
            
            expected_entropy = 1.5
            actual_entropy = df.iloc[0]["ShannonEntropy_NoGaps"]
            assert abs(actual_entropy - expected_entropy) < 0.001
            
        finally:
            alignment_file.unlink()
            if result_file and result_file.exists():
                result_file.unlink()
    
    def test_compute_conservation_scores_file_not_exists(self):
        """Test error handling for non-existent files."""
        nonexistent_file = Path("/tmp/nonexistent_alignment.fasta")
        result = compute_conservation_scores(nonexistent_file)
        assert result is None
    
    def test_compute_conservation_scores_empty_file(self):
        """Test error handling for empty alignment files."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
            # Create empty file
            pass
        
        try:
            result = compute_conservation_scores(Path(tmp.name))
            assert result is None
        finally:
            Path(tmp.name).unlink()
    
    def test_compute_conservation_scores_malformed_fasta(self):
        """Test error handling for malformed FASTA files."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
            tmp.write("This is not a valid FASTA file")
        
        try:
            result = compute_conservation_scores(Path(tmp.name))
            assert result is None
        finally:
            Path(tmp.name).unlink()
    
    def test_compute_conservation_scores_single_sequence(self):
        """Test behavior with single sequence (edge case)."""
        sequences = ["ACGT"]
        alignment_file = self.create_test_alignment_file(sequences)
        
        try:
            result_file = compute_conservation_scores(alignment_file)
            assert result_file is not None
            
            import pandas as pd
            df = pd.read_csv(result_file)
            
            # Single sequence should have entropy = 0 (perfectly conserved)
            assert len(df) == 4
            assert all(df["ShannonEntropy_WithGaps"] == 0.0)
            assert all(df["ShannonEntropy_NoGaps"] == 0.0)
            
        finally:
            alignment_file.unlink()
            if result_file and result_file.exists():
                result_file.unlink()
    
    def test_compute_conservation_scores_custom_output_path(self):
        """Test specifying custom output file path."""
        sequences = ["AAAA", "AAAA"]
        alignment_file = self.create_test_alignment_file(sequences)
        
        with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as tmp:
            custom_output = Path(tmp.name)
        
        try:
            result_file = compute_conservation_scores(alignment_file, custom_output)
            assert result_file == custom_output
            assert custom_output.exists()
            
            import pandas as pd
            df = pd.read_csv(custom_output)
            assert len(df) == 4
            
        finally:
            alignment_file.unlink()
            if custom_output.exists():
                custom_output.unlink()
    
    @patch('comparative_genomics_pipeline.service.biopython_service.path_config')
    def test_compute_conservation_for_all_msas(self, mock_path_config):
        """Test batch processing of multiple MSA files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            mock_path_config.MSA_OUTPUT_DIR = tmpdir_path
            mock_path_config.CONSERVATION_OUTPUT_DIR = tmpdir_path / "conservation"
            
            # Create multiple test MSA files
            sequences1 = ["AAAA", "AAAA"]
            sequences2 = ["CCCC", "CCCC"]
            
            msa1 = self.create_test_alignment_file(sequences1, "gene1.fasta")
            msa2 = self.create_test_alignment_file(sequences2, "gene2.fasta")
            
            # Move files to mock directory
            (tmpdir_path / "gene1.fasta").write_text(msa1.read_text())
            (tmpdir_path / "gene2.fasta").write_text(msa2.read_text())
            
            try:
                result = compute_conservation_for_all_msas()
                assert result is True
                
                # Check that conservation files were created
                conservation_dir = tmpdir_path / "conservation"
                assert conservation_dir.exists()
                
                csv_files = list(conservation_dir.glob("*.csv"))
                assert len(csv_files) == 2
                
            finally:
                msa1.unlink()
                msa2.unlink()
    
    def test_shannon_entropy_edge_cases(self):
        """Test Shannon entropy calculation edge cases."""
        # Test with gaps only
        sequences = ["----", "----", "----"]
        alignment_file = self.create_test_alignment_file(sequences)
        
        try:
            result_file = compute_conservation_scores(alignment_file)
            assert result_file is not None
            
            import pandas as pd
            df = pd.read_csv(result_file)
            
            # All gaps: with_gaps should be 0, no_gaps should handle empty case
            assert len(df) == 4
            assert all(df["ShannonEntropy_WithGaps"] == 0.0)  # All gaps = perfect conservation
            
        finally:
            alignment_file.unlink()
            if result_file and result_file.exists():
                result_file.unlink()