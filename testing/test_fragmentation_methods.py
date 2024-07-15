"""
Test file for the fragmentation_methods module.
"""
import pytest
from alphafragment.classes import Domain, Protein
from alphafragment.fragmentation_methods import validate_fragmentation_parameters, merge_overlapping_domains, check_valid_cutpoint, recursive_fragmentation

@pytest.mark.parametrize(
    "domains, expected",
    [
        # No Overlap
        ([Domain(1, 1, 5, 'TYPE'), Domain(2, 6, 10, 'TYPE')],
         [Domain(1, 1, 5, 'TYPE'), Domain(2, 6, 10, 'TYPE')]),
        # Partial Overlap
        ([Domain(1, 1, 5, 'TYPE'), Domain(2, 4, 10, 'TYPE')],
         [Domain(1, 1, 10, 'TYPE')]),
        # Complete Overlap
        ([Domain(1, 1, 10, 'TYPE'), Domain(2, 2, 9, 'TYPE')],
         [Domain(1, 1, 10, 'TYPE')]),
        # Multiple Merging
        ([Domain(1, 1, 5, 'TYPE'), Domain(2, 4, 8, 'TYPE'), Domain(3, 7, 10, 'TYPE')],
         [Domain(1, 1, 10, 'TYPE')]),
        # Touching but not overlapping
        ([Domain(1, 1, 5, 'TYPE'), Domain(2, 6, 10, 'TYPE')],
         [Domain(1, 1, 5, 'TYPE'), Domain(2, 6, 10, 'TYPE')]),
        # Order Independence
        ([Domain(2, 4, 8, 'TYPE'), Domain(1, 1, 5, 'TYPE'), Domain(3, 7, 10, 'TYPE')],
         [Domain(1, 1, 10, 'TYPE')]),
    ]
)
def test_merge_overlapping_domains(domains, expected):
    """
    Tests the merge_overlapping_domains function with various arrangements of domain
    overlaps to ensure domains are merged correctly or left separate as appropriate.
    """
    result = merge_overlapping_domains(domains)
    assert len(result) == len(expected), "Number of resulting domains does not match expected"
    for res_domain, exp_domain in zip(result, expected):
        assert res_domain.id == exp_domain.id, f"Domain id mismatch: expected {exp_domain.id}, got {res_domain.id}"
        assert res_domain.start == exp_domain.start, f"Start mismatch: expected {exp_domain.start}, got {res_domain.start}"
        assert res_domain.end == exp_domain.end, f"End mismatch: expected {exp_domain.end}, got {res_domain.end}"
        assert res_domain.type == exp_domain.type, f"Type mismatch: expected {exp_domain.type}, got {res_domain.type}"

@pytest.mark.parametrize(
    "res, domains, sequence_end, expected",
    [
        # Within Domain - not allowed
        (7, [Domain(1, 5, 10, 'TYPE')], 15, False),
        # At domain end, allowed
        (11, [Domain(1, 5, 10, 'TYPE')], 15, True),
        # At domain start, allowed
        (5, [Domain(1, 5, 10, 'TYPE')], 15, True),
        # Between Domains
        (11, [Domain(1, 1, 10, 'TYPE'), Domain(2, 11, 20, 'TYPE')], 20, True),
        # Cut before sequence start
        (-1, [], 10, False),
        # Cut at sequence start
        (0, [], 10, True),
        # Cut beyond sequence end
        (12, [], 10, False),
        # Cut at sequence end
        (11, [], 10, True),
        # End of one domain is start of another (= 1 res overlap, so not allowed)
        (10, [Domain(1, 5, 10, 'TYPE'), Domain(2, 10, 15, 'TYPE')], 15, False),
    ]
)
def test_check_valid_cutpoint(res, domains, sequence_end, expected):
    """
    Tests the check_valid_cutpoint function with various residue positions and
    domain configurations, verifying that the function correctly determines if a
    cut is valid.
    """
    result = check_valid_cutpoint(res, domains, sequence_end)
    assert result == expected, f"Expected {expected} but got {result} for residue {res} with sequence end {sequence_end}"

@pytest.mark.parametrize(
    "domains, length, overlap",
    [
        # Basic Functionality: No domains, simple split
        ([], {'min': 10, 'ideal': 15, 'max': 20}, {'min': 1, 'ideal': 2, 'max': 3}),
        # Domain Handling: No split through domain
        ([Domain(1, 5, 10, 'TYPE')], {'min': 5, 'ideal': 7, 'max': 10}, {'min': 1, 'ideal': 2, 'max': 3}),
        # Overlap Rules: Ideal overlap not possible
        ([], {'min': 5, 'ideal': 7, 'max': 10}, {'min': 1, 'ideal': 3, 'max': 4}),
    ]
)
def test_recursive_fragmentation_lengths_and_overlaps(domains, length, overlap):
    """
    Tests the recursive_fragmentation function to ensure that all fragments meet
    the length and overlap requirements specified by the parameters.
    """
    result = recursive_fragmentation(Protein("Protein1", "example_acc_id", 'A'*20),
                                     domains, 0, length, overlap)

    # Check if any fragments were generated
    for i in range(len(result)):
        # Check fragment length
        start, end = result[i]
        fragment_length = end - start
        assert length['min'] <= fragment_length <= length['max'], f"Fragment length {fragment_length} out of bounds ({length['min']}, {length['max']}) at index {i}"

        # Check overlap with next fragment if not the last one
        if i < len(result) - 1:
            next_start = result[i + 1][0]
            actual_overlap = start + fragment_length - next_start
            assert overlap['min'] <= actual_overlap <= overlap['max'], f"Overlap {actual_overlap} out of bounds ({overlap['min']}, {overlap['max']}) between fragments {i} and {i+1}"

@pytest.mark.parametrize("length, overlap, expected_error",
    [
        # max_len < min_len
        ({'min': 60, 'ideal': 55, 'max': 50}, {'min': 1, 'ideal': 2, 'max': 5}, ValueError),
        # max overlap < min_overlap
        (None, {'min': 2, 'ideal': 2, 'max': 0}, ValueError),
        # max overlap = min_len
        ({'min': 10, 'ideal': 15, 'max': 20}, {'min': 1, 'ideal': 2, 'max': 10}, ValueError),
        # ideal overlap not between min and max overlap
        (None, {'min': 1, 'ideal': 7, 'max': 5}, ValueError),
        # length values not integers
        ({'min': 10.5, 'ideal': 15, 'max': 20}, None, TypeError),
        # overlap values not integers
        (None, {'min': 1.5, 'ideal': 2, 'max': 5}, TypeError),
        # length values string
        ({'min': '10', 'ideal': '15', 'max': '20'}, None, TypeError),
        # overlap values string
        (None, {'min': '1', 'ideal': '2', 'max': '5'}, TypeError),
        # missing keys in length
        ({'min': 10, 'max': 20}, None, ValueError),
        # missing keys in overlap
        (None, {'min': 1, 'ideal': 2}, ValueError),
        # length not a dictionary
        (20, None, TypeError),
        # overlap not a dictionary
        (None, 5, TypeError),
        # min_len < 0
        ({'min': -10, 'ideal': 15, 'max': 20}, None, ValueError),
        # min_len == 0
        ({'min': 0, 'ideal': 15, 'max': 20}, None, ValueError),
        # min_overlap < 0
        (None, {'min': -1, 'ideal': 2, 'max': 5}, ValueError)
    ]
)
def test_validation_error_handling(length, overlap, expected_error):
    """
    Tests that the validate_fragmentation_parameters function raises correct
    errors for invalid length and overlap parameters.
    """
    if length is None:
        length = {'min': 10, 'ideal': 15, 'max': 20}
    if overlap is None:
        overlap = {'min': 1, 'ideal': 2, 'max': 5}
    
    try:
        validate_fragmentation_parameters(Protein("Protein1", "example_acc_id", 'A'*20),
                                          length, overlap)
    except expected_error:
        pass
    else:
        pytest.fail(f"Expected {expected_error} but no error was raised for input {length}, {overlap}")

@pytest.mark.parametrize("length, overlap",
    [
        # ideal_len == min_len
        ({'min': 10, 'ideal': 10, 'max': 20}, {'min': 1, 'ideal': 2, 'max': 5}),
        # ideal_len == max_len
        ({'min': 10, 'ideal': 20, 'max': 20}, {'min': 1, 'ideal': 2, 'max': 5}),
        # ideal_overlap == min_overlap
        ({'min': 10, 'ideal': 15, 'max': 20}, {'min': 1, 'ideal': 1, 'max': 5}),
        # ideal_overlap == max_overlap
        ({'min': 10, 'ideal': 15, 'max': 20}, {'min': 1, 'ideal': 5, 'max': 5}),
    ]
)
def test_valid_parameters_no_error(length, overlap):
    """
    Tests that the validate_fragmentation_parameters function does not raise
    an error for valid length and overlap parameters where the ideal values
    are equal to the min or max values.
    """
    try:
        validate_fragmentation_parameters(Protein("Protein1", "example_acc_id", 'A'*20),
                                          length, overlap)
    except Exception as e:
        pytest.fail(f"Unexpected error {e} raised for input {length}, {overlap}")
