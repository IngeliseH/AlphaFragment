"""
Test file for the fragment_file_creation module.
"""
from unittest.mock import mock_open, patch
import pytest
from alphafragment.classes import Protein
from alphafragment.fragment_file_creation import get_protein_combinations, output_fastas, output_pulldown

# Define sample proteins and list as a fixture
protein1 = Protein(name="ProteinA", accession_id="O123", sequence="ABCDEFGHIJKLMNOP", fragment_list=[(0,7), (6,16)])
protein2 = Protein(name="ProteinB", accession_id="P456", sequence="QRSTUVWXYZ", fragment_list=[(0,2), (1,10)])
@pytest.fixture
def proteins():
    """
    Fixture to return a list of sample Protein objects.
    """
    return [protein1, protein2]

@pytest.mark.parametrize("method, one_protein, combinations_csv, expected_output, expected_error", [
    # Valid all v all method
    ("all", None, None,
     [(protein1, protein1), (protein1, protein2), (protein2, protein2)],
     None),
    # Valid one v all method
    ("one", "ProteinA", None,
     [(protein1, protein1), (protein1, protein2)],
     None),
    # Valid specific combinations method
    ("specific", None, "dummy/path.csv",
     [(protein1, protein2)],
     None),
    # Specific combinations method with missing CSV
    ("specific", None, None, None, ValueError),
    # One v all method with missing target protein
    ("one", None, None, None, ValueError),
    # Invalid target protein name
    ("one", "ProteinC", None, None, ValueError),
    # Invalid method
    ("invalid", None, None, None, ValueError),
])
def test_get_protein_combinations(proteins, method, one_protein, combinations_csv, expected_output, expected_error):
    """
    Test the get_protein_combinations function with different methods and inputs.
    """
    if expected_error:
        with pytest.raises(expected_error):
            get_protein_combinations(proteins, method, combinations_csv=combinations_csv, one_protein=one_protein)
    else:
        mock_csv_data = "ProteinA,ProteinB\n"
        with patch("alphafragment.fragment_file_creation.open", mock_open(read_data=mock_csv_data)), \
            patch("csv.reader", return_value=[['ProteinA', 'ProteinB']]):
            combinations = get_protein_combinations(proteins, method, combinations_csv=combinations_csv, one_protein=one_protein)
            assert len(combinations) == len(expected_output)
            assert combinations == expected_output, f"Expected: {expected_output}, got: {combinations}"

@pytest.mark.parametrize("method, combinations_csv, one_protein, expected_calls", [
    # All v all method - expecting 10 files (3 self pairs for each, 4 ProteinA-ProteinB pairs)
    ("all", None, None, 10),
    # One v all method - expecting 7 files (3 self pairs for ProteinA, 4 ProteinA-ProteinB pairs)
    ("one", None, "ProteinA", 7),
    # Specific combinations method - specifying ProteinA-ProteinB so expecting 4 files (4 ProteinA-ProteinB pairs)
    ("specific", "dummy/path.csv", None, 4),
])
def test_output_fastas_creates_files(proteins, method, combinations_csv, one_protein, expected_calls):
    """
    Test the output_fastas function to ensure the correct number of files are created, with different methods.
    """
    mock_csv_data = "ProteinA,ProteinB\n"
    with patch("alphafragment.fragment_file_creation.os.makedirs"), \
         patch("alphafragment.fragment_file_creation.open", mock_open(read_data=mock_csv_data)) as mocked_file:
        output_fastas(proteins, None, method, combinations_csv=combinations_csv, one_protein=one_protein)
        handle = mocked_file()
        handle.write.assert_called(), "Expected write calls, but none were made"
        assert handle.write.call_count == expected_calls, f"Expected {expected_calls} file write operation, got {handle.write.call_count}"

def test_output_pulldown_creates_correct_files(proteins):
    """
    Test the output_pulldown function to ensure both a FASTA file and the text file (with the 'all' combination method)
    are correctly created.
    """
    with patch("alphafragment.fragment_file_creation.open", mock_open()) as mocked_file:
        output_pulldown(proteins, method="all")
        handle = mocked_file()

        # Check the number of file creation (open) calls
        open_calls = mocked_file.call_args_list
        assert len(open_calls) >= 2, f"Expected at least 2 files to be created, but found {len(open_calls)}: {open_calls}"
        
        # Extract file names from the open calls and ensure they are the correct files
        text_file_call = open_calls[0][0][0]  # First argument of the first open call (FASTA file)
        fasta_file_call = open_calls[1][0][0]   # First argument of the second open call (Text file)
        assert text_file_call.endswith('.txt'), f"Expected second file to be a text file, got: {text_file_call}"
        assert fasta_file_call.endswith('.fasta'), f"Expected first file to be a FASTA file, got: {fasta_file_call}"
        
        # Check the number of write calls
        handle.write.assert_called()
        write_calls = handle.write.call_args_list
        assert len(write_calls) >= 2, f"Expected at least 2 write calls, but got {len(write_calls)}"
        
        # Validate the content that would be written to the text file
        text_content = write_calls[0][0][0]  # Content written to the text file
        assert "O123,1-7;O123,1-7" in text_content, f"Expected 'O123,1-7;O123,1-7' in text output, got: {text_content}"
        # Ensure the correct number of lines are written in the text file
        text_lines = text_content.split('\n')
        assert len(text_lines) == 10, f"Expected 10 lines of output in the text file, got {len(text_lines)}"

        # Validate the content that would be written to the FASTA file
        fasta_content = write_calls[1][0][0]  # Content written to the FASTA file
        assert fasta_content.startswith(">"), f"Expected FASTA header, but got: {fasta_content}"

        