import pytest
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

from comparative_genomics_pipeline.service.biopython_service import compute_conservation_scores


class TestMSAProcessing:
    """Test suite for Multiple Sequence Alignment processing functions."""
    
    def create_test_alignment_file(self, sequences, descriptions=None, filename="test_alignment.fasta"):
        """Helper method to create test MSA files with optional descriptions."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
            for i, seq in enumerate(sequences):
                desc = descriptions[i] if descriptions else f"seq_{i}"
                tmp.write(f">{desc}\n{seq}\n")
            return Path(tmp.name)
    
    def create_alignment_object(self, sequences, descriptions=None):
        """Helper to create Bio.Align.MultipleSeqAlignment object."""
        records = []
        for i, seq in enumerate(sequences):
            desc = descriptions[i] if descriptions else f"seq_{i}"
            record = SeqRecord(Seq(seq), id=desc, description="")
            records.append(record)
        return MultipleSeqAlignment(records)
    
    def test_msa_file_validation_valid_file(self):
        """Test MSA file validation with valid FASTA file."""
        sequences = ["ACGT", "ACGT", "ACGT"]
        alignment_file = self.create_test_alignment_file(sequences)
        
        try:
            result = compute_conservation_scores(alignment_file)
            assert result is not None
            assert result.exists()
        finally:
            alignment_file.unlink()
            if result and result.exists():
                result.unlink()
    
    def test_msa_file_validation_file_not_exists(self):
        """Test error handling for non-existent MSA files."""
        nonexistent_file = Path("/tmp/nonexistent_msa.fasta")
        result = compute_conservation_scores(nonexistent_file)
        assert result is None
    
    def test_msa_file_validation_not_a_file(self):
        """Test error handling when path is not a file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            directory_path = Path(tmpdir)
            result = compute_conservation_scores(directory_path)
            assert result is None
    
    def test_msa_parsing_valid_fasta(self):
        """Test parsing of valid FASTA alignment file."""
        sequences = ["ATCG", "ATCG", "TTCG"]
        descriptions = ["human", "chimp", "mouse"]
        alignment_file = self.create_test_alignment_file(sequences, descriptions)
        
        try:
            result = compute_conservation_scores(alignment_file)
            assert result is not None
            
            # Verify CSV output contains expected data
            import pandas as pd
            df = pd.read_csv(result)
            assert len(df) == 4  # 4 positions
            assert "Position" in df.columns
            assert "ShannonEntropy_WithGaps" in df.columns
            assert "ShannonEntropy_NoGaps" in df.columns
            
        finally:
            alignment_file.unlink()
            if result and result.exists():
                result.unlink()
    
    def test_msa_parsing_malformed_fasta(self):
        """Test error handling for malformed FASTA files."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
            tmp.write("This is not a valid FASTA file\n")
            tmp.write("No proper headers or sequences\n")
        
        try:
            result = compute_conservation_scores(Path(tmp.name))
            assert result is None
        finally:
            Path(tmp.name).unlink()
    
    def test_msa_parsing_empty_fasta(self):
        """Test error handling for empty FASTA files."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
            pass  # Create empty file
        
        try:
            result = compute_conservation_scores(Path(tmp.name))
            assert result is None
        finally:
            Path(tmp.name).unlink()
    
    def test_msa_parsing_single_sequence(self):
        """Test parsing alignment with single sequence."""
        sequences = ["ATCGATCG"]
        alignment_file = self.create_test_alignment_file(sequences)
        
        try:
            result = compute_conservation_scores(alignment_file)
            assert result is not None
            
            import pandas as pd
            df = pd.read_csv(result)
            assert len(df) == 8  # 8 positions
            # Single sequence should have zero entropy (perfectly conserved)
            assert all(df["ShannonEntropy_WithGaps"] == 0.0)
            assert all(df["ShannonEntropy_NoGaps"] == 0.0)
            
        finally:
            alignment_file.unlink()
            if result and result.exists():
                result.unlink()
    
    def test_msa_alignment_length_validation(self):
        """Test validation of alignment length."""
        # Test with zero-length alignment (edge case)
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp:
            tmp.write(">seq1\n\n")  # Empty sequence
            tmp.write(">seq2\n\n")  # Empty sequence
        
        try:
            result = compute_conservation_scores(Path(tmp.name))
            assert result is None  # Should fail due to zero length
        finally:
            Path(tmp.name).unlink()
    
    def test_msa_alignment_length_consistency(self):
        """Test that all sequences have same length in alignment."""
        # BioPython AlignIO should handle length validation, but test our usage
        sequences = ["ATCG", "ATCG", "ATCG"]  # All same length
        alignment_file = self.create_test_alignment_file(sequences)
        
        try:
            result = compute_conservation_scores(alignment_file)
            assert result is not None
            
            import pandas as pd
            df = pd.read_csv(result)
            assert len(df) == 4  # Should process all 4 positions
            
        finally:
            alignment_file.unlink()
            if result and result.exists():
                result.unlink()
    
    def test_msa_column_extraction(self):
        """Test extraction of individual alignment columns."""
        # Test specific column patterns
        sequences = [
            "AAAA",  # All A's
            "CCCC",  # All C's
            "GGGG",  # All G's
            "TTTT"   # All T's
        ]
        alignment_file = self.create_test_alignment_file(sequences)
        
        try:
            result = compute_conservation_scores(alignment_file)
            assert result is not None
            
            import pandas as pd
            df = pd.read_csv(result)
            
            # Each position should have maximum entropy (4 different nucleotides)
            expected_entropy = 2.0  # log2(4) = 2.0
            for _, row in df.iterrows():
                assert abs(row["ShannonEntropy_NoGaps"] - expected_entropy) < 0.001
                
        finally:
            alignment_file.unlink()
            if result and result.exists():
                result.unlink()
    
    def test_msa_gap_handling(self):
        """Test proper handling of gaps in alignment columns."""
        sequences = [
            "A-CG",  # Gap in position 2
            "ATCG",  # No gaps
            "A-CG",  # Gap in position 2
            "ATCG"   # No gaps
        ]
        alignment_file = self.create_test_alignment_file(sequences)
        
        try:
            result = compute_conservation_scores(alignment_file)
            assert result is not None
            
            import pandas as pd
            df = pd.read_csv(result)
            
            # Position 2 should have different entropy with/without gaps
            pos2 = df[df["Position"] == 2].iloc[0]
            # With gaps: A(0), T(2), -(2) = mixed
            # Without gaps: T(2) only = perfect conservation
            assert pos2["ShannonEntropy_WithGaps"] > pos2["ShannonEntropy_NoGaps"]
            assert pos2["ShannonEntropy_NoGaps"] == 0.0  # Perfect conservation without gaps
            
        finally:
            alignment_file.unlink()
            if result and result.exists():
                result.unlink()
    
    def test_msa_species_names_parsing(self):
        """Test parsing of species names from FASTA headers."""
        sequences = ["ATCG", "ATCG", "ATCG"]
        species_names = ["Homo_sapiens", "Pan_troglodytes", "Mus_musculus"]
        alignment_file = self.create_test_alignment_file(sequences, species_names)
        
        try:
            result = compute_conservation_scores(alignment_file)
            assert result is not None
            
            # Function should process regardless of species names
            import pandas as pd
            df = pd.read_csv(result)
            assert len(df) == 4
            
        finally:
            alignment_file.unlink()
            if result and result.exists():
                result.unlink()
    
    def test_msa_protein_sequences(self):
        """Test processing of protein sequence alignments."""
        # Test with amino acid sequences - need identical sequences for zero entropy
        protein_sequences = [
            "MKLLVVS",
            "MKLLVVS", 
            "MKLLVVT"  # One difference at position 7
        ]
        alignment_file = self.create_test_alignment_file(protein_sequences)
        
        try:
            result = compute_conservation_scores(alignment_file)
            assert result is not None
            
            import pandas as pd
            df = pd.read_csv(result)
            assert len(df) == 7  # 7 amino acid positions
            
            # Position 7 should have some entropy (S vs T)
            pos7 = df[df["Position"] == 7].iloc[0]
            assert pos7["ShannonEntropy_NoGaps"] > 0
            
            # Positions 1-6 should be perfectly conserved (all identical)
            for pos in [1, 2, 3, 4, 5, 6]:
                pos_data = df[df["Position"] == pos].iloc[0]
                assert pos_data["ShannonEntropy_NoGaps"] == 0.0
                
        finally:
            alignment_file.unlink()
            if result and result.exists():
                result.unlink()
    
    def test_msa_large_alignment(self):
        """Test processing of larger alignments."""
        # Create alignment with more sequences and longer length
        base_sequence = "ATCGATCGATCGATCGATCG"  # 20 positions
        sequences = []
        
        # Create 10 sequences with occasional differences
        for i in range(10):
            if i == 0:
                sequences.append(base_sequence)
            else:
                # Introduce one change per sequence at different positions
                seq_list = list(base_sequence)
                change_pos = i % len(seq_list)
                seq_list[change_pos] = 'N'  # Change to N
                sequences.append(''.join(seq_list))
        
        alignment_file = self.create_test_alignment_file(sequences)
        
        try:
            result = compute_conservation_scores(alignment_file)
            assert result is not None
            
            import pandas as pd
            df = pd.read_csv(result)
            assert len(df) == 20  # 20 positions
            
            # Should successfully process large alignment
            assert df is not None
            assert len(df) > 0
            
        finally:
            alignment_file.unlink()
            if result and result.exists():
                result.unlink()
    
    @patch('comparative_genomics_pipeline.service.biopython_service.AlignIO.read')
    def test_msa_alignio_error_handling(self, mock_alignio_read):
        """Test error handling for AlignIO.read failures."""
        sequences = ["ATCG", "ATCG"]
        alignment_file = self.create_test_alignment_file(sequences)
        
        # Mock AlignIO.read to raise an exception
        mock_alignio_read.side_effect = ValueError("Mock AlignIO error")
        
        try:
            result = compute_conservation_scores(alignment_file)
            assert result is None  # Should handle AlignIO errors gracefully
            
        finally:
            alignment_file.unlink()
    
    def test_msa_utf8_encoding(self):
        """Test handling of UTF-8 encoded FASTA files."""
        sequences = ["ATCG", "ATCG", "ATCG"]
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta', encoding='utf-8') as tmp:
            for i, seq in enumerate(sequences):
                tmp.write(f">sequence_{i}\n{seq}\n")
            alignment_file = Path(tmp.name)
        
        try:
            result = compute_conservation_scores(alignment_file)
            assert result is not None
            
        finally:
            alignment_file.unlink()
            if result and result.exists():
                result.unlink()