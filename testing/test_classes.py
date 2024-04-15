"""
This test file contains tests for the classes in the AFprep_func package.
"""
import pytest
from AFprep_func.classes import Domain, Protein, ProteinSubsection  # Adjust the import as necessary.

# Tests for Domain
def test_domain_initialization():
    """
    Test that the Domain class initializes
    """
    domain = Domain("D1", 10, 20, "TYPE")
    assert domain.num == "D1"
    assert domain.start == 10
    assert domain.end == 20
    assert domain.type == "TYPE"

    # Test that using the Domain class with start/end below 0 raises a ValueError
    with pytest.raises(ValueError):
        Domain("D2", -1, 20, "TYPE")
    # Test that using the Domain class with end before start raises a ValueError
    with pytest.raises(ValueError):
        Domain("D3", 30, 20, "TYPE")

def test_domain_str():
    """
    Test that the __str__ method of the Domain class works as expected
    """
    domain = Domain("D1", 10, 20, "HELIX")
    assert str(domain) == "HELIXD1 (10, 20)"

def test_domain_equality():
    """
    Test that the __eq__ method of the Domain class works as expected
    """
    domain1 = Domain("D1", 10, 20, "TYPE")
    domain2 = Domain("D1", 10, 20, "TYPE")
    domain3 = Domain("D1", 10, 21, "TYPE")
    assert domain1 == domain2
    assert domain1 != domain3

def test_domain_same_start_end():
    """
    Test that using the Domain class with start == end raises a ValueError
    """
    with pytest.raises(ValueError):
        Domain("D1", 5, 5, "TYPE")

# Tests for Protein
def test_protein_initialization():
    """
    Test that the Protein class initializes
    """
    protein = Protein("Protein1", "P12345", "MKKLLPT", 0, 7)
    assert protein.name == "Protein1"
    assert protein.accession_id == "P12345"
    assert protein.sequence == "MKKLLPT"
    assert protein.first_res == 0
    assert protein.last_res == 7

def test_protein_add_domain():
    """
    Test that the add_domain method of the Protein class works as expected
    """
    protein = Protein("Protein1", "P12345", "MKKLLPT")
    domain = Domain("D1", 2, 5, "TYPE")
    protein.add_domain(domain)
    assert protein.domain_list == [domain]

    # Test that using the add_domain function with a non Domain class object raises a ValueError
    with pytest.raises(ValueError):
        protein.add_domain("not a domain")

def test_protein_empty_sequence():
    """
    Test Protein class initialization with an empty sequence
    """
    protein = Protein("TestProtein", "P00001", "")
    assert protein.last_res == 0

def test_protein_non_sequential_fragments():
    """
    Test that the add_fragment method of the Protein class raises a ValueError
    when non-sequential fragments are added
    """
    protein = Protein("TestProtein", "P00001", "MKKLLPT")
    protein.add_fragment(1, 3)
    with pytest.raises(ValueError):
        protein.add_fragment(0, 2)  # Tests adding an earlier fragment after a later one

# Tests for ProteinSubsection
def test_protein_subsection_initialization():
    """
    Test that the ProteinSubsection class initializes
    """
    parent_protein = Protein("Protein1", "P12345", "MKKLLPT")
    subsection = ProteinSubsection(parent_protein, 0, 3)
    assert subsection.name == parent_protein.name
    assert subsection.sequence == "MKK"

    # Test that using the ProteinSubsection class with start > end raises a ValueError
    with pytest.raises(ValueError):
        ProteinSubsection(parent_protein, 3, 2)
    # Test that using the ProteinSubsection class with start < 0 raises a ValueError
    with pytest.raises(ValueError):
        ProteinSubsection(parent_protein, -1, 5)
    # Test that using the ProteinSubsection class with end > parent_protein.last_res raises a ValueError
    with pytest.raises(ValueError):
        ProteinSubsection(parent_protein, 0, 10)

    # Test start and end end exactly at the bounds of the parent sequence
    try:
        subsection = ProteinSubsection(parent_protein, 0, len(parent_protein.sequence))
        assert subsection.sequence == parent_protein.sequence
    except ValueError:
        pytest.fail("Unexpected ValueError for a valid end boundary.")
