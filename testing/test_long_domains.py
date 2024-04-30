"""
Test file for the long_domains module
"""
import pytest
from alphafragment.classes import Protein, Domain
from alphafragment.long_domains import handle_long_domains

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
    overlap = {'min': 3, 'ideal': 5, 'max': 10}
    subsections, fragments = handle_long_domains(protein, min_len, max_len, overlap)

    assert len(subsections) == subsection_num, f"Unexpected number of subsections, expected {subsection_num}, got {len(subsections)} with subsections {subsections} and fragments {fragments}"
    assert len(fragments) == long_domain_num, f"Unexpected number of long domain fragments, expected {long_domain_num}, got {len(fragments)} with fragments {fragments} and subsections {subsections}"

    # Check that fragments fit within expected ranges
    for fragment in fragments:
        start, end = fragment
        assert max_len <= (end - start + 1), f"Fragment length too short to contain long domain, shouldn't have been created. Expected length > {max_len}, got length {end - start + 1} for fragment {fragment} in list {fragments}"
        assert start >= 0, f"Fragment start position out of bounds, expected start >= 0, got start {start} for fragment {fragment} in list {fragments}"
        assert end > start, f"Fragment end position must be after fragment start, got start {start} and end {end} for fragment {fragment} in list {fragments}"
        assert end <= len(protein_sequence), f"Fragment end position out of bounds, expected end <= {len(protein_sequence)}, got end {end} for fragment {fragment} in list {fragments}"

    # Check that subsections fit within expected ranges
    for subsection in subsections:
        subsection_len = subsection.last_res - subsection.first_res + 1
        assert subsection_len >= min_len, f"Subsection too short, should have been merged with adjacent long domain containing fragment. Expected length >= {min_len}, got length {subsection_len} for subsection {subsection} in list {subsections}"
        assert subsection.first_res >= 0, f"Subsection start position out of bounds, expected start >= 0, got start {subsection.first_res} for subsection {subsection} in list {subsections}"
        assert subsection.last_res <= len(protein_sequence), f"Subsection end position out of bounds, expected end <= {len(protein_sequence)}, got end {subsection.last_res} for subsection {subsection} in list {subsections}"

    # Check overlap between consecutive fragments and subsections
    for subsection in subsections:
        fragments.append((subsection.first_res, subsection.last_res))
    fragments_and_subsections = sorted(fragments, key=lambda x: x[0])
    for region in range(len(fragments_and_subsections) - 1):
        start, end = fragments_and_subsections[region]
        next_start, next_end = fragments_and_subsections[region + 1]
        actual_overlap = start + (end - start + 1) - next_start
        assert overlap['min'] <= actual_overlap <= overlap['max'], f"Fragment overlap out of bounds, expected overlap between ({start}, {end}) and ({next_start}, {next_end}) to be between {overlap['min']} and {overlap['max']}, got overlap {actual_overlap} for fragments/subsections {fragments_and_subsections}"

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
    min_len = 20
    max_len = 30
    overlap = {'min': 3, 'ideal': 5, 'max': 10}
    subsections, fragments = handle_long_domains(protein, min_len, max_len, overlap)

    assert len(subsections) == 0, f"Unexpected number of subsections, expected 0, got {len(subsections)} with subsections {subsections} and fragments {fragments}"
    assert len(fragments) == 2, f"Unexpected number of long domain fragments, expected 2, got {len(fragments)} with fragments {fragments} and subsections {subsections}"

    # Check that fragments fit within expected ranges
    for fragment in fragments:
        start, end = fragment
        assert max_len <= (end - start + 1), f"Fragment {fragment} length too short to contain long domain, shouldn't have been created. Fragments = {fragments}, domains = {protein.domain_list}"
        assert start >= 0, f"Fragment start position out of bounds, expected start >= 0, got start {start}"
        assert end > start, f"Fragment end position must be after fragment start, got start {start} and end {end} for fragment {fragment} in fragments {fragments}"
        assert end <= len(protein_sequence), f"Fragment end position out of bounds, expected end <= {len(protein_sequence)}, got end {end} for fragment {fragment} in fragments {fragments}"
