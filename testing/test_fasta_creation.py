import pytest
from unittest.mock import mock_open, patch
from AFprep_func.classes import Protein
from AFprep_func.fragment_fasta_creation import read_protein_combinations, generate_protein_combination_files

@pytest.mark.parametrize("file_content,expected_output", [
    # Valid protein combinations
    ("ProteinA,ProteinB\nProteinC,ProteinD\n", [("ProteinA", "ProteinB"), ("ProteinC", "ProteinD")]),
    # Empty fields for protein 1, 2 or both - should be ignored
    ("ProteinA,\n,ProteinD\n,\n", []),
    # Empty file
    ("", []),
    # Extra columns
    ("ProteinA,ProteinB,ProteinC\n", [("ProteinA", "ProteinB")]),
])
def test_read_protein_combinations(file_content, expected_output):
    with patch("builtins.open", mock_open(read_data=file_content)):
        result = read_protein_combinations("dummy/path.csv")
        assert result == expected_output, f"Expected: {expected_output}, got: {result}"

def test_read_protein_combinations_nonexistent_file():
    with pytest.raises(FileNotFoundError):
        read_protein_combinations("nonexistent/path.csv"), "Expected FileNotFoundError for nonexistent file."

# Sample Protein objects
protein1 = Protein(name="ProteinA", accession_id="", sequence="ABCDEFGHIJKLMNOP", fragment_list=[(0,7), (6,16)])
protein2 = Protein(name="ProteinB", accession_id="", sequence="QRSTUVWXYZ", fragment_list=[(0,2), (1,10)])
proteins = [protein1, protein2]

@pytest.mark.parametrize("method,combinations_csv,one_protein,expected_error", [
    # Valid all v all method
    ("all", None, None, None),
    # Valid one v all method
    ("one", None, "ProteinA", None),
    # Valid specific combinations method
    ("specific", "dummy/path.csv", None, None),
    # Specific combinations method with missing file
    ("specific", None, None, ValueError),
    # One v all method with missing target protein
    ("one", None, None, ValueError),
    # Invalid target protein name
    ("one", None, "ProteinC", ValueError),
    # Invalid method
    ("invalid", None, None, ValueError),
])
def test_generate_protein_combination_files_input_validation(method, combinations_csv, one_protein, expected_error):
    if expected_error:
        with pytest.raises(expected_error):
            generate_protein_combination_files(proteins, method, combinations_csv, one_protein)
    else:
        with patch("os.makedirs") as mocked_makedirs, patch("builtins.open", mock_open()):
            generate_protein_combination_files(proteins, method, combinations_csv, one_protein)
            assert mocked_makedirs.called, "Expected os.makedirs to be called."

def test_generate_protein_combination_files_creates_files():
    with patch("os.makedirs") as mocked_makedirs, \
         patch("builtins.open", mock_open()) as mocked_file:
        generate_protein_combination_files(proteins, method="all")
        mocked_makedirs.assert_called(), "Expected os.makedirs to be called."
        # Check calls to open to ensure correct number of files would be created - 2 proteins each with 2 fragments, results in 12 combinations including self combinations
        assert mocked_file.call_count == 12, f"Expected 12 calls to open, got {mocked_file.call_count}, {mocked_file.call_args_list}"

# Ensure content correctness in file write operations
def test_file_content_written_correctly():
    with patch("os.makedirs"), \
         patch("builtins.open", mock_open()) as mocked_file:
        generate_protein_combination_files([protein1], method="all")
        handle = mocked_file()
        # Check that the file was written with the correct content
        handle.write.assert_any_call('>ProteinA_F1+ProteinA_F2\nABCDEFG:GHIJKLMNOP'), f"Expected '>ProteinA_F1+ProteinA_F2\nMKTIIs:TIISA\n', got {handle.write.call_args_list}"