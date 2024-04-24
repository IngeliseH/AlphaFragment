"""
Test file for classes in the AFprep_func package.
"""
import pytest
from AFprep_func.classes import Domain, Protein, ProteinSubsection  # Adjust the import as necessary.

# Tests for Domain
def test_domain_initialization():
    """
    Test that the Domain class initializes. This includes checks for:
    - Regular initialization with a start and end.
    - Initialization where the start and end are the same, representing a single residue.
    - Initialization with start/end values below 0, which should raise a ValueError.
    - Initialization where the end is before the start, which should also raise a ValueError.
    """
    domain = Domain("D1", 10, 20, "TYPE")
    assert domain.id == "D1", f"Domain id did not initialize correctly, expected 'D1', got {domain.num}"
    assert domain.start == 10, f"Domain start did not initialize correctly, expected 10, got {domain.start}"
    assert domain.end == 20, f"Domain end did not initialize correctly, expected 20, got {domain.end}"
    assert domain.type == "TYPE", f"Domain type did not initialize correctly, expected 'TYPE', got {domain.type}"
    
    # Test initialization of a domain where start and end are the same (useful
    # for marking specific points of interest to see where these fall in fragments)
    single_residue_domain = Domain("SingleResidue", 5, 5, "POINT")
    assert single_residue_domain.start == 5, f"Single residue domain start did not initialize correctly, expected 5, got {single_residue_domain.start}"
    assert single_residue_domain.end == 5, f"Single residue domain end did not initialize correctly, expected 5, got {single_residue_domain.end}"
    assert single_residue_domain.type == "POINT", f"Single residue domain type did not initialize correctly, expected 'POINT', got {single_residue_domain.type}"
    
    # Test that using the Domain class with start/end below 0 raises a ValueError
    with pytest.raises(ValueError):
        Domain("D2", -1, 20, "TYPE"), "Domain with start below 0 did not raise a ValueError"
    # Test that using the Domain class with end before start raises a ValueError
    with pytest.raises(ValueError):
        Domain("D3", 30, 20, "TYPE"), "Domain with end before start did not raise a ValueError"

def test_domain_str():
    """
    Test that the __str__ method of the Domain class works as expected
    """
    domain = Domain("D1", 10, 20, "HELIX")
    assert str(domain) == "HELIXD1 (10, 20)", f"Expected str(domain) to give 'HELIXD1 (10, 20)', got {str(domain)}"

def test_domain_equality():
    """
    Test that the __eq__ method of the Domain class works as expected
    """
    domain1 = Domain("D1", 10, 20, "TYPE")
    domain2 = Domain("D1", 10, 20, "TYPE")
    domain3 = Domain("D1", 10, 21, "TYPE")
    assert domain1 == domain2, "Equal domains did not return True for equality check"
    assert domain1 != domain3, "Different domains returned True for equality check"

# Tests for Protein
def test_protein_initialization():
    """
    Test that the Protein class initializes
    """
    protein = Protein("Protein1", "P12345", "MKKLLPT", 0, 7)
    assert protein.name == "Protein1", f"Protein name did not initialize correctly, expected 'Protein1', got {protein.name}"
    assert protein.accession_id == "P12345", f"Protein accession ID did not initialize correctly, expected 'P12345', got {protein.accession_id}"
    assert protein.sequence == "MKKLLPT", f"Protein sequence did not initialize correctly, expected 'MKKLLPT', got {protein.sequence}"
    assert protein.first_res == 0, f"Protein first_res did not initialize correctly, expected 0, got {protein.first_res}"
    assert protein.last_res == 7, f"Protein last_res did not initialize correctly, expected 7, got {protein.last_res}"

def test_protein_add_domain():
    """
    Test that the add_domain method of the Protein class works as expected
    """
    protein = Protein("Protein1", "P12345", "MKKLLPT")
    domain = Domain("D1", 2, 5, "TYPE")
    protein.add_domain(domain)
    assert protein.domain_list == [domain], f"add_domain did not add the domain to the protein's domain list, expected [Domain('D1', 2, 5, 'TYPE')], got {protein.domain_list}"

    # Test that using the add_domain function with a non Domain class object raises a ValueError
    with pytest.raises(ValueError):
        protein.add_domain("not a domain"), "add_domain did not raise a ValueError for a non-Domain object"

def test_protein_empty_sequence():
    """
    Test Protein class initialization with an empty sequence
    """
    protein = Protein("TestProtein", "P00001", "")
    assert protein.last_res == None, f"Protein last_res did not initialize correctly for an empty sequence, expected 0, got {protein.last_res}"

@pytest.mark.parametrize("fragment, expected_exception, message", [
    # Valid sequential fragments
    ((30, 40), None, ""),
    # Overlap with previous fragment
    ((15, 60), None, ""),
    # Start same as previous
    ((10, 50), None, ""),
    # Fragment start is greater than end
    ((30, 25), ValueError, "Start and end must be positive integers, and start must be less than end."),
    # Invalid input type - not a tuple (or two integers)
    (30, ValueError, "Invalid arguments. Provide (start, end) as either two arguments or a tuple."),
    # Invalid input type - not both integers
    (('a', 60), ValueError, "Start and end must be positive integers, and start must be less than end."),
    # Non-sequential fragment
    ((1, 5), ValueError, "Start of the new fragment must be greater than the start of the previous fragment.")
])
def test_add_fragment(fragment, expected_exception, message):
    protein = Protein("TestProtein", "fake_id", "A"*100)
    #Add initial (valid) fragment
    protein.add_fragment(10, 20)
    if expected_exception:
        with pytest.raises(ValueError) as exc_info:
            protein.add_fragment(fragment), "Expected ValueError not raised."
        assert message in str(exc_info.value), f"Unexpected error message: {exc_info.value}. Expected: {message}"
    else:
        protein.add_fragment(fragment)
        # Verify by checking the last added fragment
        assert protein.fragment_list[-1] == fragment, "Fragment was not added correctly."

# Tests for ProteinSubsection
def test_protein_subsection_initialization():
    """
    Test that the ProteinSubsection class initializes
    """
    parent_protein = Protein("Protein1", "P12345", "MKKLLPT")
    subsection = ProteinSubsection(parent_protein, 0, 3)
    assert subsection.name == parent_protein.name, f"ProteinSubsection name did not initialize correctly, expected 'Protein1', got {subsection.name}"
    assert subsection.sequence == "MKKL", f"ProteinSubsection sequence did not initialize correctly, expected 'MKK', got {subsection.sequence}"

    # Test that using the ProteinSubsection class with start > end raises a ValueError
    with pytest.raises(ValueError):
        ProteinSubsection(parent_protein, 3, 2), "ProteinSubsection with start > end did not raise a ValueError"
    # Test that using the ProteinSubsection class with start < 0 raises a ValueError
    with pytest.raises(ValueError):
        ProteinSubsection(parent_protein, -1, 5), "ProteinSubsection with start < 0 did not raise a ValueError"
    # Test that using the ProteinSubsection class with end > parent_protein.last_res raises a ValueError
    with pytest.raises(ValueError):
        ProteinSubsection(parent_protein, 0, 10), "ProteinSubsection with end > parent_protein.last_res did not raise a ValueError"

    # Test start and end exactly at the bounds of the parent sequence
    try:
        subsection = ProteinSubsection(parent_protein, 0, parent_protein.last_res)
        assert subsection.sequence == parent_protein.sequence, f"ProteinSubsection sequence did not initialize correctly at the bounds of the parent sequence, expected {parent_protein.sequence}, got {subsection.sequence}"
    except ValueError:
        pytest.fail("Unexpected ValueError for a valid end boundary")
