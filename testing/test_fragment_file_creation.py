"""
Test file for the fragment_file_creation module.
"""

import pytest
from unittest.mock import mock_open, patch
from AFprep_func.classes import Protein
from AFprep_func.fragment_file_creation import get_protein_combinations, output_fastas, output_pulldown

# Define sample proteins and list as a fixture
protein1 = Protein(name="ProteinA", accession_id="O123", sequence="ABCDEFGHIJKLMNOP", fragment_list=[(0,7), (6,16)])
protein2 = Protein(name="ProteinB", accession_id="P456", sequence="QRSTUVWXYZ", fragment_list=[(0,2), (1,10)])
@pytest.fixture
def proteins():
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
        with patch("AFprep_func.fragment_file_creation.open", mock_open(read_data=mock_csv_data)), \
            patch("csv.reader", return_value=[['ProteinA', 'ProteinB']]):
            combinations = get_protein_combinations(proteins, method, combinations_csv=combinations_csv, one_protein=one_protein)
            assert len(combinations) == len(expected_output)
            assert combinations == expected_output, f"Expected: {expected_output}, got: {result}"

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
    with patch("AFprep_func.fragment_file_creation.os.makedirs"), \
         patch("AFprep_func.fragment_file_creation.open", mock_open(read_data=mock_csv_data)) as mocked_file:
        output_fastas(proteins, None, method, combinations_csv=combinations_csv, one_protein=one_protein)
        handle = mocked_file()
        handle.write.assert_called(), "Expected write calls, but none were made"
        assert handle.write.call_count == expected_calls, f"Expected {expected_calls} file write operation, got {handle.write.call_count}"

def test_output_pulldown_creates_correct_file(proteins):
    """
    Test the output_pulldown function to ensure the correct file is created (with the 'all' combination method).
    """
    with patch("AFprep_func.fragment_file_creation.open", mock_open()) as mocked_file:
        output_pulldown(proteins, method="all")
        handle = mocked_file()
        handle.write.assert_called()
        calls = handle.write.call_args_list
        assert len(calls) > 0, "Expected write calls, but none were made"
        lines = calls[0][0][0].split('\n') # Split output into lines
        print(lines)
        first_line = lines[0]
        assert first_line == "O123,1-7;O123,1-7", f"Incorrect first line of output, got: {first_line}"
        # Check that the correct number of lines are written
        # Expecting 10 - 3 for each self pair ((1,1), (1,2), (2,2)), 4 for ProteinA-ProteinB ((1,1), (1,2), (2,1), (2,2))
        assert len(lines) == 10, f"Expected 10 lines of output, got {len(lines)}: {calls}"
