"""
This test file contains tests for the alphafold_db_domain_identification module.
"""
import pytest
import requests_mock
from requests.exceptions import Timeout, ConnectionError as RequestsConnectionError
from alphafragment.classes import Domain
from alphafragment.alphafold_db_domain_identification import read_afdb_json, find_domain_by_res, find_domains_from_pae

@pytest.fixture
def mock_requests():
    """
    Provides a requests mocking fixture for use in test functions to simulate API
    responses.
    """
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
    """
    Test the fetching of PAE data from a mocked API endpoint to verify handling
    of different responses.
    """
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
    """
    Tests the function's response to various HTTP error statuses to ensure that
    it properly handles errors like 404 and 500.
    """
    mock_url = 'https://alphafold.ebi.ac.uk/files/AF-invalid_id-F1-predicted_aligned_error_v4.json'
    mock_requests.get(mock_url, status_code=status_code)

    result = read_afdb_json("invalid_id", "v4")
    assert result is None, f"Expected None for HTTP error {status_code}, got {result}"

@pytest.mark.parametrize("exception", [Timeout, RequestsConnectionError])
def test_request_exceptions_handling(mock_requests, exception):
    """
    Tests the handling of different request exceptions such as timeouts and connection
    errors to ensure the function gracefully handles these scenarios.
    """
    mock_url = 'https://alphafold.ebi.ac.uk/files/AF-valid_id-F1-predicted_aligned_error_v4.json'
    mock_requests.get(mock_url, exc=exception)

    result = read_afdb_json("valid_id", "v4")
    assert result is None, f"Expected None for exception {exception}, got {result}"

@pytest.mark.parametrize("accession_id", [None, 'na', '', 'NA'])
def test_read_afdb_json_with_invalid_accession_id(accession_id):
    """
    Test that read_afdb_json returns None when provided with an invalid or non-applicable accession_id.
    """
    result = read_afdb_json(accession_id)
    assert result is None, f"Expected None for invalid accession_id '{accession_id}', but got {result}"

@pytest.mark.parametrize("domains, residue, expected", [
    # Test with no domains
    ([], 5, None),
    # Test with one domain, residue inside
    ([Domain("AF1", 1, 10, 'AF')], 5, Domain("AF1", 1, 10, 'AF')),
    # Test with multiple domains, residue in second domain
    ([Domain("AF1", 1, 5, 'AF'), Domain("AF2", 6, 10, 'AF')], 7, Domain("AF2", 6, 10, 'AF')),
    # Test with multiple domains, residue not in any domain
    ([Domain("AF1", 1, 5, 'AF'), Domain("AF2", 6, 10, 'AF')], 11, None)
])
def test_find_domain_by_res(domains, residue, expected):
    """
    Tests the function that determines domain membership by residue number by
    checking several scenarios.
    """
    result = find_domain_by_res(domains, residue)
    assert result == expected

@pytest.mark.parametrize("pae, method, expected_exception", [
    # PAE input that is not a matrix (list of integers)
    ([1, 2, 3, 4], 'cautious', ValueError),
    # PAE input that is not square
    ([[1, 2], [3, 4], [5, 6]], 'cautious', ValueError),
    # PAE contains non-numeric values
    ([[1, 2], [3, 'a']], 'cautious', ValueError),
    # Invalid method
    ([[1, 2], [3, 4]], 'invalid', ValueError)
])
def test_error_handling(pae, method, expected_exception):
    """
    Ensures that the domain finding function correctly raises exceptions for various general input errors, including incorrect PAE matrices and methods.
    """
    with pytest.raises(expected_exception):
        find_domains_from_pae(pae, method)

@pytest.mark.parametrize("custom_params, expected_exception", [
    # Complete and correct custom parameters
    ({"res_dist_cutoff": 5, "close_pae_val": 1, "further_pae_val": 2}, None),
    # Custom method with non-dict parameters
    ([1, 2], ValueError),
    # Custom method missing necessary params
    ({'res_dist_cutoff': 5}, ValueError),
    # Non-numeric value in parameters
    ({"res_dist_cutoff": "five", "close_pae_val": 0.5, "further_pae_val": 2}, ValueError),
])
def test_handling_custom_parameters(custom_params, expected_exception):
    """
    Tests the handling of custom parameters in the find_domains_from_pae function,
    both in correct cases and where errors should be raised.
    """
    pae = [[1, 2], [3, 4]]  # Common PAE matrix used in all cases
    method = 'custom'        # Common method used in all cases

    if expected_exception:
        with pytest.raises(expected_exception):
            find_domains_from_pae(pae, method, custom_params)
    else:
        # Test should pass without raising an exception
        result = find_domains_from_pae(pae, method, custom_params)
        assert isinstance(result, list), f"Expected a list of domains for valid inputs, got {type(result)}"


@pytest.mark.parametrize("pae, expected_length", [
    # A very small matrix
    ([[1]], 0),
    # A very large matrix
    ([[3]*5000 for _ in range(5000)], 1),
    # Matrix with only high values (no domains expected)
    ([[30]*10 for _ in range(10)], 0),
    # 4x4 matrix with low values, but within ignored distance so no domains expected
    ([[1]*4 for _ in range(4)], 0),
    # Negative values in PAE matrix (expect all residues to be grouped into one domain)
    ([[-1]*10 for _ in range(10)], 1)
])
def test_edge_cases_input_structure(pae, expected_length):
    """
    Tests the domain finding function with PAE data structures that represent edge
    cases to ensure domains are identified or not identified as expected.
    """
    domains = find_domains_from_pae(pae)
    # simplify pae if large so failure print statemtents not overwhelming
    if len(pae) > 5:
        pae = f"length {len(pae)}, starting {pae[0][:10]}..."
    assert len(domains) == expected_length, f"Expected {expected_length} domains, got {len(domains)} domains, {domains} for pae {pae}"

@pytest.mark.parametrize("pae, expected_domain_length", [
    # Whole matrix with very low values (all in one domain)
    ([[1]*10 for _ in range(10)], 10),
    # Sparse confidence test (whole matrix in one domain but only confident on edges)
    ([[0,   111, 111, 111, 111, 111, 0,   111, 111, 111],
      [111, 0,   111, 111, 111, 111, 111, 111, 111, 111],
      [111, 111, 0,   111, 111, 111, 111, 111, 111, 111],
      [111, 111, 111, 0,   111, 111, 111, 111, 111, 111],
      [111, 111, 111, 111, 0,   111, 111, 111, 111, 111],
      [111, 111, 111, 111, 111, 0,   111, 111, 111, 111],
      [0,   111, 111, 111, 111, 111, 0,   111, 111, 111],
      [111, 111, 111, 111, 111, 111, 111, 0,   111, 111],
      [111, 111, 111, 111, 111, 111, 111, 111, 0,   111],
      [111, 111, 111, 111, 111, 111, 111, 111, 111, 0]], 7),
    # Correctly extends domain if confident residue pairs overlap
    ([[0,   111, 111, 111, 111, 0,   111, 111, 111, 111, 111],
      [111, 0,   111, 111, 111, 111, 111, 0,   111, 111, 111],
      [111, 111, 0,   111, 111, 111, 111, 111, 111, 111, 111],
      [111, 111, 111, 0,   111, 111, 111, 111, 111, 111, 111],
      [111, 111, 111, 111, 0,   111, 111, 111, 111, 111, 111],
      [0,   111, 111, 111, 111, 0,   111, 111, 111, 111, 111],
      [111, 111, 111, 111, 111, 111, 0,   111, 111, 111, 111],
      [111, 0,   111, 111, 111, 111, 111, 0,   111, 111, 111],
      [111, 111, 111, 111, 111, 111, 111, 111, 0,   111, 111],
      [111, 111, 111, 111, 111, 111, 111, 111, 111, 0,   111],
      [111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 0]], 8),
    # Non-symmetric PAE matrix, still should find overlapping confident regions
    ([[0,   111, 111, 111, 111, 0,   111, 111, 111, 111],
      [111, 0,   111, 111, 111, 111, 111, 111, 111, 111],
      [111, 111, 0,   111, 111, 111, 111, 111, 111, 111],
      [111, 111, 111, 0,   111, 111, 111, 111, 111, 111],
      [111, 111, 111, 111, 0,   111, 111, 111, 111, 111],
      [111, 111, 111, 111, 111, 0,   111, 111, 111, 111],
      [111, 111, 111, 111, 111, 111, 0,   111, 111, 111],
      [111, 111, 111, 111, 111, 111, 111, 0,   111, 111],
      [111, 111, 111, 0,   111, 111, 111, 111, 0,   111],
      [111, 111, 111, 111, 111, 111, 111, 111, 111, 0]], 9),
    # Still correctly extends domain if overlap is on boundary
    ([[0,   111, 111, 111, 111, 0,   111, 111, 111, 111, 111],
      [111, 0,   111, 111, 111, 111, 111, 111, 111, 111, 111],
      [111, 111, 0,   111, 111, 111, 111, 111, 111, 111, 111],
      [111, 111, 111, 0,   111, 111, 111, 111, 111, 111, 111],
      [111, 111, 111, 111, 0,   111, 111, 111, 111, 111, 111],
      [111, 111, 111, 111, 111, 0,   111, 111, 111, 111, 111],
      [111, 111, 111, 111, 111, 111, 0,   111, 111, 111, 111],
      [111, 111, 111, 111, 111, 111, 111, 0,   111, 111, 111],
      [111, 111, 111, 111, 111, 111, 111, 111, 0,   111, 111],
      [111, 111, 111, 111, 111, 111, 111, 111, 111, 0,   111],
      [111, 111, 111, 111, 111, 0,   111, 111, 111, 111, 0]], 11),
])
def test_correct_domain_finding(pae, expected_domain_length):
    """
    Verifies that domains are correctly identified and extended.
    """
    domains = find_domains_from_pae(pae)
    assert len(domains) == 1, f"Expected 1 domain, got {len(domains)} domains, {domains}"
    assert domains[0].end - domains[0].start + 1 == expected_domain_length, f"Expected domain length {expected_domain_length}, got {domains[0].end - domains[0].start + 1}"

def test_different_methods():
    """
    Compares the results of domain identification using different methods ('cautious'
    vs 'definite') to ensure that each method yields expected results/differences.
    """
    pae_matrix = [[0,  1,  1,  1,  1,  1,  2,  11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                  [1,  0,  1,  1,  1,  1,  2,  11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                  [1,  1,  0,  1,  1,  1,  2,  11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                  [1,  1,  1,  0,  1,  1,  2,  11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                  [1,  1,  1,  1,  0,  1,  2,  11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                  [1,  1,  1,  1,  1,  0,  2,  11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                  [2,  2,  2,  2,  2,  2,  0,  11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                  [11, 11, 11, 11, 11, 11, 11, 0,  11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                  [11, 11, 11, 11, 11, 11, 11, 11, 0,  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                  [12, 12, 12, 12, 12, 12, 12, 12, 12, 0,  1,  1,  1,  1,  1,  5,  5,  5,  5,  5,  5],
                  [12, 12, 12, 12, 12, 12, 12, 12, 12, 1,  0,  1,  1,  1,  1,  5,  5,  5,  5,  5,  5],
                  [12, 12, 12, 12, 12, 12, 12, 12, 12, 1,  1,  0,  1,  1,  1,  5,  5,  5,  5,  5,  5],
                  [12, 12, 12, 12, 12, 12, 12, 12, 12, 1,  1,  1,  0,  1,  1,  5,  5,  5,  5,  5,  5],
                  [12, 12, 12, 12, 12, 12, 12, 12, 12, 1,  1,  1,  1,  0,  1,  5,  5,  5,  5,  5,  5],
                  [12, 12, 12, 12, 12, 12, 12, 12, 12, 1,  1,  1,  1,  1,  0,  5,  5,  5,  5,  5,  5],
                  [12, 12, 12, 12, 12, 12, 12, 12, 12, 5,  5,  5,  5,  5,  5,  0,  1,  1,  1,  1,  1],
                  [12, 12, 12, 12, 12, 12, 12, 12, 12, 5,  5,  5,  5,  5,  5,  1,  0,  1,  1,  1,  1],
                  [12, 12, 12, 12, 12, 12, 12, 12, 12, 5,  5,  5,  5,  5,  5,  1,  1,  0,  1,  1,  1],
                  [12, 12, 12, 12, 12, 12, 12, 12, 12, 5,  5,  5,  5,  5,  5,  1,  1,  1,  0,  1,  1],
                  [12, 12, 12, 12, 12, 12, 12, 12, 12, 5,  5,  5,  5,  5,  5,  1,  1,  1,  1,  0,  1],
                  [12, 12, 12, 12, 12, 12, 12, 12, 12, 5,  5,  5,  5,  5,  5,  1,  1,  1,  1,  1,  0]
                  ]
    cautious_domains = find_domains_from_pae(pae_matrix, 'cautious')
    definite_domains = find_domains_from_pae(pae_matrix, 'definite')
    assert len(definite_domains) > len(cautious_domains), f"Expected more domains with 'definite' method, got {len(cautious_domains)} domains with 'cautious' method and {len(definite_domains)} domains with 'definite' method"
