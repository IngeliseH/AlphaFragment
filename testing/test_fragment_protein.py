"""
Test file for the fragment_protein module.
"""
import pytest
from functions.fragment_protein import fragment_protein
from functions.classes import Protein, Domain

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
    overlap = {'min': 5, 'ideal': 10, 'max': 15}
    fragments = fragment_protein(protein, min_len, max_len, overlap)
    for i in range(len(fragments) - 1):
        # Check correct fragment length
        fragment_length = fragments[i][1] - fragments[i][0]
        assert min_len <= fragment_length <= max_len, f"Fragment length of fragment {fragments[i]} out of bounds, for fragment index {i} in fragment list {fragments}"
        # Check correct overlap between fragments
        actual_overlap = fragments[i][1] - fragments[i + 1][0]
        assert overlap['min'] <= actual_overlap <= overlap['max'], f"Overlap of {actual_overlap} out of bounds between fragments {i} and {i+1} in fragment list {fragments}"

def test_long_domain_handling():
    """
    Test that the function correctly handles long domains
    """
    protein = Protein(name="Protein1", accession_id="example_acc_id",
                      sequence='A'*275, domain_list=[Domain(1, 20, 150, 'TYPE')])
    overlap = {'min': 5, 'ideal': 7, 'max': 10}
    fragments = fragment_protein(protein, 50, 100, overlap)
    # Check correct number of fragments
    assert len(fragments) == 3, f"Expected 3 fragments, got {len(fragments)}"
    # Check correct overlap between fragments
    for i in range(len(fragments) - 1):
        actual_overlap = fragments[i][1] - fragments[i + 1][0]
        assert overlap['min'] <= actual_overlap <= overlap['max'], f"Overlap of {actual_overlap} out of bounds between fragments {i} and {i+1} in fragment list {fragments}"

def test_very_small_protein():
    """
    Test that the function correctly handles very small proteins (shorter than min fragment length
    """
    protein = Protein(name="Protein1", accession_id="example_acc_id",
                      sequence='A'*10, domain_list=[])
    fragments = fragment_protein(protein, min_len=20, max_len=50, overlap={'min':0, 'ideal':5, 'max':5})
    assert len(fragments) == 1, f"Expected single fragment for very small protein, got {len(fragments)} fragments"
    assert fragments[0] == (0, 10), f"Expected fragment to cover full length of protein (0, 10) for very small protein, got {fragments[0]}"

@pytest.mark.parametrize("min_len, max_len, min_overlap, ideal_overlap, max_overlap", [
    # min_len > max_len
    (60, 50, 0, 5, 5),
    # overlap less than min_overlap
    (20, 50, 20, 15, 25),
    # overlap greater than max_overlap
    (20, 50, 10, 40, 30),
    # max_overlap is greater than min_len
    (20, 50, 5, 10, 25)
])
def test_invalid_parameters(min_len, max_len, min_overlap, ideal_overlap, max_overlap):
    """
    Test that the function raises the correct exception for invalid parameters
    """
    overlap = {'min': min_overlap, 'ideal': ideal_overlap, 'max': max_overlap}
    protein = Protein(name="Protein1", accession_id="example_acc_id",
                      sequence='A'*100, domain_list=[Domain(1, 10, 90, 'TYPE')])
    with pytest.raises(ValueError):
        fragment_protein(protein, min_len, max_len, overlap), f"Expected ValueError for input {min_len, max_len, overlap}, got {fragment_protein(protein, min_len, max_len, overlap)}"

def test_invalid_protein_input():
    """
    Test that the function raises the correct exception for invalid protein input
    """
    with pytest.raises(TypeError):
        fragment_protein("not a protein"), "Expected TypeError for invalid protein input"
