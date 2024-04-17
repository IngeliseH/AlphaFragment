"""
Test file for the fragment_protein module.
"""
import pytest
from AFprep_func.fragment_protein import fragment_protein
from AFprep_func.classes import Protein, Domain

@pytest.mark.parametrize(
    "protein", [
        # Basic sequence fragmentation
        (Protein(name="BasicProtein", accession_id="basic_acc", sequence='A'*100,
                 domain_list=[Domain(1, 0, 4, 'TYPE'), Domain(2, 5, 9, 'TYPE'), Domain(3, 10, 14, 'TYPE'), Domain(4, 15, 19, 'TYPE'), Domain(5, 20, 24, 'TYPE')])),
        # Very long sequence with no domains
        (Protein(name="LongNoDomainProtein", accession_id="long_no_domain_acc", sequence='A'*1000, domain_list=[])),
        # Sequence with many short domains
        (Protein(name="ManyShortDomainsProtein", accession_id="many_short_domains_acc", sequence='A'*100,
                 domain_list=[Domain(1, 1, 10, 'TYPE'), Domain(2, 17, 20, 'TYPE'), Domain(3, 33, 40, 'TYPE'), Domain(4, 60, 70, 'TYPE'), Domain(5, 85, 95, 'TYPE')])),
        # Very long sequence with many short domains
        (Protein(name="LongManyShortDomainsProtein", accession_id="long_many_short_domains_acc", sequence='A'*1000,
                 domain_list=[Domain(i, i*50, i*50+20, 'TYPE') for i in range(20)]))
])
def test_basic_fragmentation(protein):
    """
    Test basic functionality with various inputs
    """
    min_len = 20
    max_len = 50
    ideal_overlap = 10
    min_overlap = 5
    max_overlap = 15
    fragments = fragment_protein(protein, min_len, max_len, ideal_overlap, min_overlap, max_overlap)
    for i in range(len(fragments) - 1):
        # Check correct fragment length
        fragment_length = fragments[i][1] - fragments[i][0]
        assert min_len <= fragment_length <= max_len, f"Fragment length of fragment {fragments[i]} out of bounds, for fragment index {i} in fragment list {fragments}"
        # Check correct overlap between fragments
        actual_overlap = fragments[i][1] - fragments[i + 1][0]
        assert min_overlap <= actual_overlap <= max_overlap, f"Overlap of {actual_overlap} out of bounds between fragments {i} and {i+1} in fragment list {fragments}"

def test_long_domain_handling():
    """
    Test that the function correctly handles long domains
    """
    protein = Protein(name="Protein1", accession_id="example_acc_id",
                      sequence='A'*275, domain_list=[Domain(1, 20, 150, 'TYPE')])
    ideal_overlap = 7
    min_overlap = 5
    max_overlap = 10
    fragments = fragment_protein(protein, 50, 100, ideal_overlap, min_overlap, max_overlap)
    # Check correct number of fragments
    assert len(fragments) == 3, f"Expected 3 fragments, got {len(fragments)}"
    # Check correct overlap between fragments
    for i in range(len(fragments) - 1):
        actual_overlap = fragments[i][1] - fragments[i + 1][0]
        assert min_overlap <= actual_overlap <= max_overlap, f"Overlap of {actual_overlap} out of bounds between fragments {i} and {i+1} in fragment list {fragments}"

def test_very_small_protein():
    """
    Test that the function correctly handles very small proteins (shorter than min fragment length
    """
    protein = Protein(name="Protein1", accession_id="example_acc_id",
                      sequence='A'*10, domain_list=[])
    fragments = fragment_protein(protein, 20, 50, 5, 0, 5)
    assert len(fragments) == 1, f"Expected single fragment for very small protein, got {len(fragments)} fragments"
    assert fragments[0] == (0, 10), f"Expected fragment to cover full length of protein (0, 10) for very small protein, got {fragments[0]}"

@pytest.mark.parametrize("params,expected_exception", [
    # min_len > max_len
    ((60, 50, 5, 0, 5), ValueError),
    # overlap less than min_overlap
    ((20, 50, 15, 20, 25), ValueError),
    # overlap greater than max_overlap
    ((20, 50, 40, 10, 30), ValueError),
    # max_overlap is greater than min_len
    ((20, 50, 10, 5, 25), ValueError),
    # invalid protein input
    ("not a protein", TypeError)
])
def test_invalid_input(params, expected_exception):
    if isinstance(params, tuple):
        min_len, max_len, overlap, min_overlap, max_overlap = params
        protein = Protein(name="Protein1", accession_id="example_acc_id",
                          sequence='A'*100, domain_list=[Domain(1, 10, 90, 'TYPE')])
        with pytest.raises(expected_exception):
            fragment_protein(protein, min_len, max_len, overlap, min_overlap, max_overlap), f"Expected {expected_exception} for input {params}"
    else:
        with pytest.raises(expected_exception):
            fragment_protein(params, 20, 50, 5, 0, 5), f"Expected {expected_exception} for input {params}"