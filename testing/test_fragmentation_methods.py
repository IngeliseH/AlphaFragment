"""
Test file for the fragmentation_methods module.
"""
import pytest
from AFprep_func.classes import Domain, Protein
from AFprep_func.fragmentation_methods import merge_overlapping_domains, check_valid_cutpoint, recursive_fragmentation

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
        assert res_domain.num == exp_domain.num, f"Domain num mismatch: expected {exp_domain.num}, got {res_domain.num}"
        assert res_domain.start == exp_domain.start, f"Start mismatch: expected {exp_domain.start}, got {res_domain.start}"
        assert res_domain.end == exp_domain.end, f"End mismatch: expected {exp_domain.end}, got {res_domain.end}"
        assert res_domain.type == exp_domain.type, f"Type mismatch: expected {exp_domain.type}, got {res_domain.type}"

@pytest.mark.parametrize(
    "res, domains, sequence_end, expected",
    [
        # Within Domain - not allowed
        (5, [Domain(1, 1, 10, 'TYPE')], 10, False),
        # At domain end, allowed
        (10, [Domain(1, 1, 10, 'TYPE')], 10, True),
        # At domain start, allowed
        (1, [Domain(1, 1, 10, 'TYPE')], 10, True),
        # Between Domains
        (11, [Domain(1, 1, 10, 'TYPE'), Domain(2, 12, 20, 'TYPE')], 20, True),
        # Cut before sequence start
        (-1, [], 10, False),
        # Cut beyond sequence end
        (11, [], 10, False),
        # Cut at sequence end
        (10, [], 10, True),
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
    "domains, min_len, max_len, ideal_overlap, min_overlap, max_overlap",
    [
        # Basic Functionality: No domains, simple split
        ([], 10, 15, 2, 1, 3),
        # Domain Handling: No split through domain
        ([Domain(1, 5, 10, 'TYPE')], 5, 10, 2, 1, 3),
        # Overlap Rules: Ideal overlap not possible
        ([], 5, 10, 3, 1, 4),
    ]
)
def test_recursive_fragmentation_lengths_and_overlaps(domains, min_len, max_len,
                                                      ideal_overlap, min_overlap, max_overlap):
    """
    Tests the recursive_fragmentation function to ensure that all fragments meet
    the length and overlap requirements specified by the parameters.
    """
    result = recursive_fragmentation(Protein("Protein1", "example_acc_id", 'A'*20), domains, 0, min_len,
                                     max_len, ideal_overlap, min_overlap, max_overlap)

    # Check if any fragments were generated
    for i in range(len(result)):
        # Check fragment length
        start, end = result[i]
        fragment_length = end - start
        assert min_len <= fragment_length <= max_len, f"Fragment length {fragment_length} out of bounds ({min_len}, {max_len}) at index {i}"

        # Check overlap with next fragment if not the last one
        if i < len(result) - 1:
            next_start = result[i + 1][0]
            actual_overlap = start + fragment_length - next_start
            assert min_overlap <= actual_overlap <= max_overlap, f"Overlap {actual_overlap} out of bounds ({min_overlap}, {max_overlap}) between fragments {i} and {i+1}"

@pytest.mark.parametrize("ideal_overlap, min_overlap, max_overlap",
    [
        # max overlap < min_overlap
        (2, 2, 0),
        # max overlap = min_len
        (2, 1, 5),
        # ideal overlap not between min and max overlap
        (7, 1, 5),
])
def test_recursive_fragmentation_overlap_error(ideal_overlap, min_overlap, max_overlap):
    """
    Tests that the recursive_fragmentation function raises a ValueError when
    overlap and length constraints are not met.
    """
    with pytest.raises(ValueError):
        recursive_fragmentation(Protein("Protein1", "example_acc_id", 'A'*20),
                                [], 0, 5, 15, ideal_overlap,
                                min_overlap, max_overlap), f"Function did not raise ValueError for invalid overlap parameters. Input parameters: {ideal_overlap}, {min_overlap}, {max_overlap}"
