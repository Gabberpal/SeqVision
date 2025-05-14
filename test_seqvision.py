import unittest
from unittest.mock import patch
import os
from io import StringIO
from bio_files_processor import convert_multiline_fasta_to_oneline, parse_blast_output
from SeqVision import DNASequence, RNASequence, AminoAcidSequence, filter_fastq


class TestDNAOperations(unittest.TestCase):
    """Tests for DNA sequence operations."""

    def test_dna_complement(self):
        dna = DNASequence("ATGC")
        result = dna.complement()
        self.assertEqual(str(result), "TACG")

    def test_dna_reverse_complement(self):
        dna = DNASequence("ATGC")
        result = dna.reverse_complement()
        self.assertEqual(str(result), "GCAT")

    def test_invalid_dna_sequence(self):
        with self.assertRaises(ValueError):
            DNASequence("ATGCX")  # Invalid character 'X'


class TestRNAOperations(unittest.TestCase):
    """Tests for RNA sequence operations."""

    def test_rna_complement(self):
        rna = RNASequence("AUGC")
        result = rna.complement()
        self.assertEqual(str(result), "UACG")

    def test_rna_reverse_complement(self):
        rna = RNASequence("AUGC")
        result = rna.reverse_complement()
        self.assertEqual(str(result), "GCAU")


class TestAminoAcidOperations(unittest.TestCase):
    """Tests for protein sequence operations."""

    def test_protein_molecular_weight(self):
        protein = AminoAcidSequence("MKWV")
        result = protein.get_molecular_weight()
        self.assertEqual(result, 616.78)  # Expected molecular weight for "MKWV"


class TestFileOperations(unittest.TestCase):
    """Tests for file handling operations."""

    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    def test_convert_fasta_to_oneline(self, mock_file):
        # Simulate the FASTA conversion process
        input_fasta = "data/example_multiline_fasta.fasta"
        output_fasta = "output_fasta.fasta"

        # Mock the conversion function to simulate reading and writing files
        convert_multiline_fasta_to_oneline(input_fasta, output_fasta)

        # Check if the open function was called with the correct arguments
        mock_file.assert_called_with("filtered/output_fasta.fasta", "a")


class TestErrorHandling(unittest.TestCase):
    """Tests for error handling and invalid input."""

    def test_invalid_file_format(self):
        with self.assertRaises(FileNotFoundError):
            convert_multiline_fasta_to_oneline(
                "data/invalid_file.fasta", "filtered/output.fasta"
            )

    def test_invalid_sequence_type(self):
        with self.assertRaises(ValueError):
            DNASequence("ATGCX")  # Invalid character 'X'


if __name__ == "__main__":
    unittest.main()
