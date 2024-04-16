"""
Test file for the long_domains module
"""
import pytest
from AFprep_func.classes import Protein, Domain
from AFprep_func.long_domains import handle_long_domains

@pytest.mark.parametrize(
    "protein_domains, protein_sequence, long_domain_num, subsection_num",
    [   # Test case with one long domain at the start of the sequence
        ([Domain(1, 1, 30, 'TYPE')], 'A'*80, 1, 1),
        # Test case with one long domain in the middle of the sequence - ends are
        # not long enough to be subsections so should be merged
        ([Domain(1, 10, 50, 'TYPE')], 'A'*40, 1, 0),
        # Test case with one long domain in the middle of the sequence - ends are
        # long enough to be subsections so shouldn't be merged
        ([Domain(1, 35, 70, 'TYPE')], 'A'*110, 1, 2),
        # No domains, entire sequence is one subsection
        ([], 'A'*40, 0, 1),
        # No long domains, entire sequence is one subsection
        ([Domain(1, 5, 10, 'TYPE')], 'A'*40, 0, 1),
        # Test case with one long and one short domain close to each other and to boundaries
        # expecting merging
        ([Domain(1, 1, 10, 'TYPE'), Domain(2, 15, 45, 'TYPE')], 'A'*50, 1, 0),
        # Test case with a long domain created from two overlapping short domains
        ([Domain(1, 1, 20, 'TYPE'), Domain(2, 15, 31, 'TYPE')], 'A'*70, 1, 1)
])
def test_handle_long_domains(protein_domains, protein_sequence, long_domain_num, subsection_num):
    """
    Test for handling of long domains within a protein sequence.
    """
    protein = Protein(name='protein1', accession_id='example_acc_id',
                      sequence=protein_sequence, domain_list=protein_domains,
                      first_res=0, last_res=len(protein_sequence) - 1)
    min_len = 20
    max_len = 30
    ideal_overlap = 5
    min_overlap = 3
    max_overlap = 10
    subsections, fragments = handle_long_domains(protein, min_len, max_len,
                                                 ideal_overlap, min_overlap, max_overlap)

    assert len(subsections) == subsection_num, "Unexpected number of subsections"
    assert len(fragments) == long_domain_num, "Unexpected number of long domain fragments"

    # Check that fragments fit within expected ranges
    for fragment in fragments:
        start, end = fragment
        assert max_len <= (end - start + 1), "Fragment length too short to contain long domain, shouldn't have been created"
        assert start >= 0, "Fragment start position out of bounds"
        assert end > start, "Fragment end position must be after fragment start"
        assert end <= len(protein_sequence), "Fragment end position out of bounds"

    # Check that subsections fit within expected ranges
    for subsection in subsections:
        subsection_len = subsection.last_res - subsection.first_res + 1
        assert subsection_len >= min_len, "Subsection too short, should have been merged with adjacent long domain containing fragment"
        assert subsection.first_res >= 0, "Subsection start position out of bounds"
        assert subsection.last_res <= len(protein_sequence), "Subsection end position out of bounds"

    # Check overlap between consecutive fragments and subsections
    for subsection in subsections:
        fragments.append((subsection.first_res, subsection.last_res))
    fragments_and_subsections = sorted(fragments, key=lambda x: x[0])
    for region in range(len(fragments_and_subsections) - 1):
        start, end = fragments_and_subsections[region]
        next_start, next_end = fragments_and_subsections[region + 1]
        overlap = start + (end - start + 1) - next_start
        assert min_overlap <= overlap <= max_overlap, "Fragment overlap out of bounds"

def test_adjacent_long_domains():
    """
    Test for handling adjacent long domains that should be treated as two
    subsections (despite not meeting minimum overlap requirements.)
    """
    protein_sequence = 'A'*70
    protein = Protein(name='protein1', accession_id='example_acc_id',
                      sequence=protein_sequence,
                      domain_list=[Domain(1, 0, 31, 'TYPE'), Domain(2, 32, 70, 'TYPE')],
                      first_res=0, last_res=len(protein_sequence) - 1)
    max_len = 30
    subsections, fragments = handle_long_domains(protein, min_len=20, max_len=30,
                                                 overlap=5, min_overlap=3, max_overlap=10)

    assert len(subsections) == 0, "Unexpected number of subsections"
    assert len(fragments) == 2, "Unexpected number of long domain fragments"

    # Check that fragments fit within expected ranges
    for fragment in fragments:
        start, end = fragment
        assert max_len <= (end - start + 1), f"Fragment {fragment} length too short to contain long domain {fragments}, shouldn't have been created"
        assert start >= 0, "Fragment start position out of bounds"
        assert end > start, "Fragment end position must be after fragment start"
        assert end <= len(protein_sequence), "Fragment end position out of bounds"
