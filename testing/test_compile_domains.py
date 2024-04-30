"""
Test file for the compile_domains module
"""
from unittest.mock import patch
import pytest
import pandas as pd
from functions.classes import Protein, Domain
from functions.compile_domains import compile_domains

valid_protein_data = pd.DataFrame({
    'name': ['TestProtein'],
    'domains': [[(60, 70), (80, 90), (95, 100)]]})
invalid_protein_data = pd.DataFrame({'name': ['TestProtein']})

@pytest.mark.parametrize("uniprot,alphafold,manual,protein_data,expected_count,expected_print", [
    # Uniprot only
    (True, False, False, None, 1, ""),
    # AlphaFold only
    (False, True, False, None, 2, "2 domains found in AlphaFold structure for TestProtein: [Domain(id=domain2, start=20, end=30, domain_type='AF'), Domain(id=domain3, start=40, end=50, domain_type='AF')]"),
    # Manual only (with valid data)
    (False, False, True, valid_protein_data, 3, ""),
    # All three (with valid data)
    (True, True, True, valid_protein_data, 6, "2 domains found in AlphaFold structure for TestProtein: [Domain(id=domain2, start=20, end=30, domain_type='AF'), Domain(id=domain3, start=40, end=50, domain_type='AF')]"),
    # Manual only, with no data provided
    (False, False, True, None, 0, "Cannot add manually specified domains - the 'df' argument is required."),
    # Manual only, with invalid data
    (False, False, True, invalid_protein_data, 0, "Cannot add manually specified domains - the 'domains' column does not exist in the dataframe."),
    # None, with valid data to check data not added when manual False
    (False, False, False, valid_protein_data, 0, ""),
    # Check that other data is still added when manual is True but data is not provided
    (True, True, True, None, 3, "Cannot add manually specified domains - the 'df' argument is required."),
])
def test_compile_domains(uniprot, alphafold, manual, protein_data, expected_count, expected_print, capsys):
    """
    Test for the compile_domains function with different combinations of data sources.
    """
    uniprot_data = [Domain('domain1', 0, 10, 'UniProt')]
    alphafold_data = [Domain('domain2', 20, 30, 'AF'), Domain('domain3', 40, 50, 'AF')]
    manual_data = [Domain('domain4', 60, 70, 'manual'), Domain('domain5', 80, 90, 'manual'), Domain('domain6', 95, 100, 'manual')]
    with patch('AFprep_func.compile_domains.find_uniprot_domains', return_value=uniprot_data), \
         patch('AFprep_func.compile_domains.read_afdb_json', return_value='dummy_pae'), \
         patch('AFprep_func.compile_domains.find_domains_from_pae', return_value=alphafold_data), \
         patch('AFprep_func.compile_domains.find_user_specified_domains', return_value=manual_data):

        protein = Protein("TestProtein", "mock_accid", "sequence")
        domains = compile_domains(protein, uniprot=uniprot, alphafold=alphafold, manual=manual, protein_data=protein_data)

        assert len(domains) == expected_count, f"Expected {expected_count} domains, got {len(domains)}, domains: {domains}"
        if expected_print:
            captured = capsys.readouterr()
            assert expected_print in captured.out, f"Expected print statement not found: expected {expected_print}, got {captured.out}"

# Test that the function raises a TypeError if the protein argument is not an instance of the Protein class
def test_compile_domains_wrong_protein():
    """
    Test that the compile_domains function raises a TypeError if the 'protein' argument is not an instance of the Protein class.
    """
    with pytest.raises(TypeError, match="The 'protein' argument must be an instance of the Protein class."):
        compile_domains("not a protein", manual = False)
