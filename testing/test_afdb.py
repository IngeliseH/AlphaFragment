import pytest
import requests_mock
from AFprep_func.alphafold_db_domain_identification import read_afdb_json
from requests.exceptions import Timeout, ConnectionError

@pytest.fixture
def mock_requests():
    with requests_mock.Mocker() as m:
        yield m

@pytest.mark.parametrize("response, expected", [
    # Empty dictionary in a list
    ([{}], None),
    # Invalid JSON string
    ("not a json", None),
    # Correct structure but missing the 'predicted_aligned_error' key
    ([{"missing_key": "no PAE data"}], None),
    # Valid response with 'predicted_aligned_error' data
    ([{"predicted_aligned_error": [[1, 2], [3, 4]]}], [[1, 2], [3, 4]]),
])
def test_fetch_pae_data(mock_requests, response, expected):
    mock_url = 'https://alphafold.ebi.ac.uk/files/AF-valid_id-F1-predicted_aligned_error_v4.json'
    
    # If the response is a string, simulate a non-JSON response by setting text, else set json
    if isinstance(response, str):
        mock_requests.get(mock_url, text=response)
    else:
        mock_requests.get(mock_url, json=response)
    
    result = read_afdb_json("valid_id", "v4")
    assert result == expected, f"Expected {expected} for response {response}, got {result}"

@pytest.mark.parametrize("status_code", [404, 500])
def test_http_error_handling(mock_requests, status_code):
    mock_url = 'https://alphafold.ebi.ac.uk/files/AF-invalid_id-F1-predicted_aligned_error_v4.json'
    mock_requests.get(mock_url, status_code=status_code)
    
    result = read_afdb_json("invalid_id", "v4")
    assert result is None, f"Expected None for HTTP error {status_code}, got {result}"

@pytest.mark.parametrize("exception", [Timeout, ConnectionError])
def test_request_exceptions_handling(mock_requests, exception):
    mock_url = 'https://alphafold.ebi.ac.uk/files/AF-valid_id-F1-predicted_aligned_error_v4.json'
    mock_requests.get(mock_url, exc=exception)
    
    result = read_afdb_json("valid_id", "v4")
    assert result is None, f"Expected None for exception {exception}, got {result}"