"""
Test file for the uniprot_fetch module
"""
from unittest.mock import patch
import pytest
import requests_mock
from requests.exceptions import ConnectionError as RequestsConnectionError
from functions.uniprot_fetch import fetch_uniprot_info, find_uniprot_domains
from functions.classes import Protein, Domain

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
     {'features': [{'type': 'DOMAIN', 'begin': '1', 'end': '50', 'description': 'Domain_1_description'},
                   {'type': 'CHAIN', 'begin': '1', 'end': '100'},  # Should be ignored
                   {'type': 'DOMAIN', 'begin': '51', 'end': '100', 'description': 'Domain_2_description'}]},
     [Domain('Domain_1_description', 0, 49, 'UniProt'), Domain('Domain_2_description', 50, 99, 'UniProt')],
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
        assert result == expected_domains, f"Failed to find expected domains {expected_domains} in test protein. Should have returned {expected_output}. Actual result was {result}"
