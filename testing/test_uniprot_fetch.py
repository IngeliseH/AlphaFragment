"""
Test file contains test cases for the uniprot_fetch module
"""
import pytest
import requests_mock
from requests.exceptions import ConnectionError as RequestsConnectionError
from unittest.mock import patch
from AFprep_func.uniprot_fetch import fetch_uniprot_info, find_uniprot_domains
from AFprep_func.classes import Protein, Domain

@pytest.mark.parametrize("accession_id, mock_url, response, status_code, expected, exception", [
    # Valid ID with empty features
    ("P12345", "https://www.ebi.ac.uk/proteins/api/features/P12345", {'features': []}, 200, {'features': []}, None),
    # Invalid ID
    ("INVALID", "https://www.ebi.ac.uk/proteins/api/features/INVALID", None, 404, None, None),
    # Server error
    ("P12345", "https://www.ebi.ac.uk/proteins/api/features/P12345", None, 500, None, None),
    # Connection error
    ("P12345", "https://www.ebi.ac.uk/proteins/api/features/P12345", None, None, None, RequestsConnectionError),
])
def test_fetch_uniprot_info(accession_id, mock_url, response, status_code, expected, exception):
    """
    Test the 'fetch_uniprot_info' function with different inputs and requests errors.
    """
    with requests_mock.Mocker() as m:
        if exception:
            m.get(mock_url, exc=exception)
        else:
            m.get(mock_url, json=response, status_code=status_code)

        # Assuming 'fetch_uniprot_info' is the function being tested from your module
        result = fetch_uniprot_info(accession_id)
        if exception:
            assert result is None, f"Expected None when an exception occurs, but got {result}"
        else:
            assert result == expected, f"Expected {expected}, but got {result}"

@pytest.mark.parametrize("accession_id, uniprot_data, expected_domains, expected_output", [
    # No data available
    ("P12345", None, None, None),
    # Valid but empty features list
    ("P12345", {'features': []}, [], []),
    # Valid data with multiple domains
    ("P12345",
     {'features': [{'type': 'DOMAIN', 'begin': '1', 'end': '50', 'description': 'Test domain 1'},
                   {'type': 'CHAIN', 'begin': '1', 'end': '100'},  # Should be ignored
                   {'type': 'DOMAIN', 'begin': '51', 'end': '100', 'description': 'Test domain 2'}]},
     [Domain('Test domain 1', 1, 50, 'DOMAIN'), Domain('Test domain 2', 51, 100, 'DOMAIN')],
     "Expected extracted domains")
])
def test_find_uniprot_domains(accession_id, uniprot_data, expected_domains, expected_output):
    """
    Test the 'find_uniprot_domains' function with different inputs.
    """
    protein = Protein(name="Test Protein", accession_id=accession_id, sequence="")
    # Mock the fetch function to return our test data without making an actual API call
    with patch('AFprep_func.uniprot_fetch.fetch_uniprot_info', return_value=uniprot_data):
        result = find_uniprot_domains(protein)
        assert result == expected_domains, f"Failed {expected_output} for {accession_id}"