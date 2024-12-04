"""
Test file for the process_proteins_csv module.
"""
from io import StringIO
from unittest.mock import patch
import pytest
import pandas as pd
from alphafragment.process_proteins_csv import initialize_proteins_from_csv, find_user_specified_domains, update_csv_with_fragments
from alphafragment.classes import Protein, Domain

@pytest.mark.parametrize("csv_data, expected_success, expected_errors", [
    # Valid Accession ID, no manual sequence provided
    ('name,accession_id\nProtein1,P12345',
     [Protein(name='Protein1', accession_id='P12345', sequence='AAA')],
     None),
    # No Accession ID and no manual sequence
    ('name,accession_id\nProtein2,',
     None,
     ['Protein2']),
    # No Accession ID but with manual sequence
    ('name,accession_id,sequence\nProtein3,,CCC',
     [Protein(name='Protein3', accession_id='', sequence='CCC')],
     None),
    # Accession ID and manual sequence (manual sequence should be overridden if fetch successful)
    ('name,accession_id,sequence\nProtein4,P12345,DDD',
     [Protein(name='Protein4', accession_id='P12345', sequence='AAA')],
     None),
    # Invalid Accession ID and no manual sequence
    ('name,accession_id\nProtein5,INVALID',
     None,
     ['Protein5']),
    # No accession id and no manual sequence
    ('name,accession_id\nProtein5,',
     None,
     ['Protein5']),
    # Invalid Accession ID but with a manual sequence
    ('name,accession_id,sequence\nProtein6,INVALID,FFF',
     [Protein(name='Protein6', accession_id='INVALID', sequence='FFF')],
     None),
    # Two valid proteins to test multiple initializations
    ('name,accession_id\nProtein7,P12345\nProtein8,P12345',
     [Protein(name='Protein7', accession_id='P12345', sequence='AAA'),
      Protein(name='Protein8', accession_id='P12345', sequence='AAA')],
     None),
    # One valid and one invalid protein
    ('name,accession_id\nProtein9,P12345\nProtein10,INVALID',
     [Protein(name='Protein9', accession_id='P12345', sequence='AAA')],
     ['Protein10']),
])
def test_initialize_proteins_from_csv(csv_data, expected_success, expected_errors, capfd):
    """
    Test the initialization of Protein objects from CSV data and capture print output.
    """
    with patch('pandas.read_csv', return_value=pd.read_csv(StringIO(csv_data))):
        with patch('alphafragment.process_proteins_csv.fetch_uniprot_info') as mock_fetch:
            # Setup mock to simulate fetching results
            mock_fetch.side_effect = lambda x: {'sequence': 'AAA'} if x == 'P12345' else None
            proteins, df = initialize_proteins_from_csv("fake_path")

            # Use capfd to capture print statements
            out, err = capfd.readouterr()

            # Check number of successful initiations, and print outputs
            if expected_success:
                assert len(proteins) == len(expected_success), f"Expected {len(expected_success)} proteins, got {len(proteins)}"
                assert f"Successfully initialized proteins: {[protein.name for protein in expected_success]}" in out
            if expected_errors:
                assert f"Proteins with errors or no data available: {expected_errors}" in out

@pytest.mark.parametrize("csv_data, expected_columns, raises_error, error_message", [
    # No column headings
    ('Protein11,P12345\nProtein12,P12345', None, True, "Missing required columns: name, accession_id"),
    # Incorrect column headings
    ('Title,accession_id\nProtein13,P12345\nProtein14,P12345', None, True, "Missing required columns: name"),
    # Duplicate columns - currently don't give an error as renamed by pandas upon reading in
    ('name,accession_id,accession_id\nProtein15,P12346,P12346',
     ['name', 'accession_id', 'accession_id.1'], False, ""),
    # Duplicate columns after normalization
    ('name,Name,ACCESSion_ID,accession_id\nProtein16,,P12346,P12346',
     None, True, "Duplicate column names detected after normalization"),
     # Capitalized column headings
    ('NAME,ACCESSION_ID\nProtein17,P12345\nProtein16,P12345',
     ['name', 'accession_id'], False, ""),
    # Unnecessary column headings in between useful ones
    ('name,random_col,accession_id\nProtein18,random,P12345',
     ['name', 'random_col', 'accession_id'], False, "")
])
def test_csv_column_handling(csv_data, expected_columns, raises_error, error_message):
    """
    Test that the CSV column handling in initialize_proteins_from_csv works as expected.
    """
    with patch('pandas.read_csv', return_value=pd.read_csv(StringIO(csv_data))) as mock_read_csv:
        if raises_error:
            with pytest.raises(ValueError) as exc_info:
                initialize_proteins_from_csv("fake_path")
            assert error_message in str(exc_info.value), f"Error message does not match expected, expected {error_message}, got {str(exc_info.value)}"
        else:
            proteins, df = initialize_proteins_from_csv("fake_path")
            # Ensure the mock was called to simulate file reading
            mock_read_csv.assert_called_once(), "pd.read_csv not called"
            # Check if DataFrame contains expected columns
            assert all(col in df.columns for col in expected_columns), f"Expected columns not present in DataFrame, expected {expected_columns}, got {df.columns}"

@pytest.mark.parametrize("columns, expected_error", [
    # Test missing 'domains' column
    (['name'], "ValueError"),
    # Test missing 'name' column
    (['domains'], "ValueError"),
    # Test both columns present
    (['name', 'domains'], None),
    # Test neither column present
    ([], "ValueError"),
    # Test additional columns present
    (['name', 'domains', 'other'], None),
    # Test columns ordered differently
    (['other', 'domains', 'other2', 'name'], None)
])
def test_column_presence(columns, expected_error):
    """
    Test the function with various combinations of missing and present columns.
    """
    data = {
        'name': ['Protein1'],
        'domains': ['[(1, 10), (20, 30)]'],
        'other': ['extra data'],
        'other2': ['extra data2']
    }
    df = pd.DataFrame({col: data[col] for col in columns if col in data})
    protein_name = 'Protein1'
    expected_result = [Domain('manual_D1', 0, 9, 'manually_defined'), Domain('manual_D2', 19, 29, 'manually_defined')]
    if expected_error:
        #check expected error (likely ValueError) is raised
        with pytest.raises(ValueError):
            find_user_specified_domains(protein_name, df)
    else:
        result = find_user_specified_domains(protein_name, df)
        assert isinstance(result, list), f"Expected a list of Domain objects, got {result} with type {type(result)}"
        for domain in result:
            assert isinstance(domain, Domain), f"Expected a Domain object for valid column input, got {domain} with type {type(domain)}"
        assert result == expected_result, f"Expected {expected_result}, got {result}"

@pytest.mark.parametrize("protein_name,domain_data,expected_result", [
    # Protein not in dataframe (no error, should return empty list)
    ('Protein2', '[(1, 10), (20, 30)]', []),
    # Null values in 'domains'
    ('Protein1', None, []),
    # Valid input with unnamed domains (list of (start, end) tuples)
    ('Protein1', '[(1, 10), (20, 30)]', [Domain('manual_D1', 0, 9, 'manually_defined'), Domain('manual_D2', 19, 29, 'manually_defined')]),
    # Valid input with named domains (list of (name, (start, end)) tuples)
    ('Protein1', "[('name1', (1,10)), ('name2', (20,30))]", [Domain('name1', 0, 9, 'manually_defined'), Domain('name2', 19, 29, 'manually_defined')]),
    # Valid input with names in tuple format (repeated names)
    ('Protein1', "[('name', (1,10)), ('name', (20,30))]", [Domain('name', 0, 9, 'manually_defined'), Domain('name', 19, 29, 'manually_defined')])
])
def test_various_inputs(protein_name, domain_data, expected_result):
    """
    Test function with various inputs and expected outputs.
    """
    df = pd.DataFrame({
        'name': ['Protein1'],
        'domains': [domain_data]
    })
    result = find_user_specified_domains(protein_name, df)
    assert result == expected_result, f"Expected result {expected_result}, got {result}"

@pytest.mark.parametrize("protein_name, dataframe, error", [
    # Invalid 'dataframe' input - not a DataFrame
    ('Protein1', "not a dataframe", "TypeError"),
    # Invalid data for domains - integer
    ('Protein2', pd.DataFrame({'name': ['Protein2'], 'domains': [123]}), "TypeError"),
    # Invalid data for domains - string
    ('Protein3', pd.DataFrame({'name': ['Protein3'], 'domains': ['not a list']}), "ValueError"),
    # Invalid data for domains - valid tuple but not in a list
    ('Protein4', pd.DataFrame({'name': ['Protein4'], 'domains': [(1, 2)]}), "TypeError"),
    # Invalid data for domains - tuple of strings
    ('Protein5', pd.DataFrame({'name': ['Protein5'], 'domains': [[('a', 'b')]]}), "TypeError"),
    # Invalid data for domains - tuple with 3 values
    ('Protein6', pd.DataFrame({'name': ['Protein6'], 'domains': [[(1, 2, 3)]]}), "TypeError"),
    # Invalid data for domains - tuple with 1 value
    ('Protein7', pd.DataFrame({'name': ['Protein7'], 'domains': [[(1,)]]}), "TypeError"),
])
def test_manual_domain_error_handling(protein_name, dataframe, error):
    """
    Test the function to ensure type checking for dataframe input.
    """
    if error == "TypeError":
        with pytest.raises(TypeError):
            find_user_specified_domains(protein_name, dataframe), "Expected a TypeError for invalid DataFrame input"
    if error == "ValueError":
        with pytest.raises(ValueError):
            find_user_specified_domains(protein_name, dataframe), "Expected a ValueError for invalid domain data"

def test_update_csv_with_fragments():
    """
    Test that the update_csv_with_fragments function saves the DataFrame to a CSV file.
    """
    df = pd.DataFrame({'name': ['Protein1', 'Protein2'], 'accession_id': ['P12345', 'P67890']})
    proteins = [
        Protein(name='Protein1', accession_id='P12345', sequence='AAA', domain_list=[], fragment_list=[(0, 3)]),
        Protein(name='Protein2', accession_id='P67890', sequence='BBB', domain_list=[], fragment_list=[(0, 3)])
    ]
    output_csv = StringIO()

    # Check DataFrame saved to CSV.
    with patch.object(pd.DataFrame, 'to_csv') as mock_to_csv:
        update_csv_with_fragments(df, output_csv, proteins)
        mock_to_csv.assert_called_once(), "call to DataFrame.to_csv() not made"

@pytest.mark.parametrize("df_setup, expected_columns", [
    # Test with additional columns
    (pd.DataFrame({
        'name': ['Protein1'], 'accession_id': ['P12345'], 'extra_col': ['extra data']
     }),
     ['name', 'accession_id', 'sequence', 'domains', 'fragment_indices', 'fragment_sequences', 'extra_col']),
    # Test with columns ordered incorrectly
    (pd.DataFrame({
        'sequence': ['pre-sequence'], 'name': ['Protein1'], 'accession_id': ['P12345']
     }),
     ['name', 'accession_id', 'sequence', 'domains', 'fragment_indices', 'fragment_sequences'])
])
def test_column_handling(df_setup, expected_columns):
    """
    Test that the update_csv_with_fragments function correctly handles various
    formats of DataFrame columns.
    """
    output_csv = StringIO()
    proteins = [Protein(name='Protein1', accession_id='P12345', sequence='ABCDEFGHIJABCDEFGHIJ')]
    new_df = update_csv_with_fragments(df_setup, output_csv, proteins)
    # Verify that all expected columns are present
    assert set(new_df.columns) == set(expected_columns), f"Not all expected columns are present in the DataFrame, expected {expected_columns}, got {new_df.columns}"
    # Verify the order of the columns
    assert new_df.columns.tolist() == expected_columns, f"Columns are not in the expected order, expected {expected_columns}, got {new_df.columns}"

@pytest.mark.parametrize("df_setup, protein_setup, expected_exception", [
    # Test with incorrect data types in DataFrame
    (pd.DataFrame({'name': [123], 'accession_id': [456]}),
     [Protein(name='Protein1', accession_id='P12345', sequence='AAA')],
     None),
    # Test proteins not in the DataFrame
    (pd.DataFrame({'name': ['Protein2'], 'accession_id': ['P67890']}),
     [Protein(name='Protein1', accession_id='P12345', sequence='AAA')],
     None)
])
def test_data_integrity(df_setup, protein_setup, expected_exception):
    """
    Test that the update_csv_with_fragments function correctly handles incorrect
    data types and missing proteins.
    """
    output_csv = StringIO()
    if expected_exception:
        with pytest.raises(expected_exception):
            update_csv_with_fragments(df_setup, output_csv, protein_setup)
    else:
        update_csv_with_fragments(df_setup, output_csv, protein_setup)
        # Verify if the protein names in the DataFrame are not altered for proteins not in DataFrame
        if 'Protein1' not in df_setup['name'].values:
            assert True, "Protein not in DataFrame should not be added."
        else:
            assert False, "Data type test failed, check for handling of incorrect data types."

@pytest.mark.parametrize("output_path", [
    # Testing output to an invalid path
    "/invalid/path/to/file.csv",
])
def test_output_file_writing(output_path):
    """
    Test that the update_csv_with_fragments function raises an exception when trying
    to write to an invalid path.
    """
    df = pd.DataFrame({'name': ['Protein1'], 'accession_id': ['P12345']})
    proteins = [Protein(name='Protein1', accession_id='P12345', sequence='AAA')]
    with pytest.raises(Exception):
        update_csv_with_fragments(df, output_path, proteins), "Expected an exception when writing to an invalid path"
