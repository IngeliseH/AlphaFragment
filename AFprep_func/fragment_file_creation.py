"""
This module contains functions for creating protein fragment pairs and files using
these, for use as input to AlphaFold or AlphaPulldown.

Functions:
    - get_protein_combinations: Generates a list of protein pairs based on the method selected.
    - output_fastas: Creates fasta files for each combination of protein fragments.
    - output_pulldown: Creates a txt file with protein fragment combinations compatible with AlphaPulldown.

Dependencies:
    - csv: For reading CSV files.
    - os: For creating directories and files.
    - combinations: For generating combinations of proteins.
"""
import csv
import os
from itertools import combinations
from AFprep_func.classes import Protein

def get_protein_combinations(proteins, method, combinations_csv, one_protein):
    """
    Generates a list of protein pairs based on the method selected.
    
    Parameters:
        - proteins (list): A list of Protein objects.
        - method (str): Which protein combinations to use - 'all' for all
          v all, 'one' for all v one, 'specific' for specific combinations
          (specified in a csv file using the combinations_csv argument).
        - combinations_csv (csv path, optional): The path to a CSV file containing
          specific protein combinations. Each row in the CSV should represent one
          combination, with two columns for the names of the two proteins, and no
          column headings. Required if method is 'specific', otherwise ignored.
        - one_protein (str, optional): Name of a protein to use for one v all
          combinations. Required if method is 'one', otherwise ignored.
    
    Returns:
        - protein_combinations (list): A list of tuples, where each tuple contains
          two Protein objects.
    
    Raises:
        - ValueError: If any of the inputs are invalid.
    """
    # Validate inputs
    if not all(isinstance(protein, Protein) for protein in proteins):
        raise ValueError("Input proteins must be Protein objects.")
    if method not in ['all', 'one', 'specific']:
        raise ValueError("Method must be 'all', 'one', or 'specific'.")
    if method == 'specific' and combinations_csv is None:
        raise ValueError("Specific combinations selected but no combinations CSV provided.")
    if method == 'one' and not one_protein:
        raise ValueError("Method 'one' selected but no target protein specified.")

    # For 'all' method, compile all combinations of proteins, including self-combinations
    elif method == 'all':
        protein_pairs = list(combinations(proteins, 2)) + [(p, p) for p in proteins]
    
    # For 'one' method, find the target protein and compile all combinations with it, including self-combinations
    elif method == 'one':
        target_protein = next((p for p in proteins if p.name == one_protein), None)
        if not target_protein:
            raise ValueError(f"Protein named {one_protein} not found among the provided proteins.")
        protein_pairs = [(target_protein, p) for p in proteins]

    # For 'specific' method, read the desired combinations from the CSV
    if method == 'specific':
        protein_dict = {protein.name: protein for protein in proteins}
        protein_pairs = []
        with open(combinations_csv, newline='') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                if row[0] and row[1]:
                # Check if both proteins specified in the row exist in the protein list
                    if row[0] in protein_dict and row[1] in protein_dict:
                        protein_pairs.append((protein_dict[row[0]], protein_dict[row[1]]))
                    else:
                        # Print error message if one or both proteins are not found
                        missing_proteins = [name for name in row[:2] if name not in protein_dict]
                        print(f"Error: Combination {row[0]}-{row[1]} not possible. Missing proteins: {', '.join(missing_proteins)}.")

    return sorted(protein_pairs, key=lambda pair: (pair[0].name, pair[1].name))

def output_fastas(proteins, method='all', combinations_csv=None, one_protein=None):
    """
    Creates a folder for each protein pair, containing fasta files for each combination
    of fragments for that pair.
    
    Parameters:
        - proteins (list): A list of Protein objects.
        - method (str, optional): Which protein combinations to use - 'all' for all
          v all, 'one' for all v one, 'specific' for specific combinations
          (specified using the combinations_csv argument). Default is 'all'.
        - combinations_csv (csv path, optional): The path to a CSV file containing
          specific protein combinations to generate. Each row in the CSV should
          represent one combination, with two columns for the names of the two proteins.
          Default is None. Required if method is 'specific', otherwise ignored.
        - one_protein (str, optional): Name of a protein to use for one v all
          combinations. Default is None. Required if method is 'one', otherwise
          ignored.
    
    Returns:
        - None, but creates folders and .fasta files for each protein pair.
    """
    protein_pairs = get_protein_combinations(proteins, method, combinations_csv, one_protein)

    for protein1, protein2 in protein_pairs:
        folder_name = f"{protein1.name}_{protein2.name}"
        os.makedirs(folder_name, exist_ok=True)
        for i, (start1, end1) in enumerate(protein1.fragment_list):
            for j, (start2, end2) in enumerate(protein2.fragment_list):
                filename = os.path.join(folder_name, f"{protein1.name}_F{i+1}+{protein2.name}_F{j+1}.fasta")
                content = f">{protein1.name}_F{i+1}+{protein2.name}_F{j+1}\n{protein1.sequence[start1:end1]}:{protein2.sequence[start2:end2]}"
                with open(filename, 'w') as file:
                    file.write(content)
                print(f"File created: {filename}")

def output_pulldown(proteins, output_name='pulldown_input.txt', method='all', combinations_csv=None, one_protein=None):
    """
    Creates a txt file with protein fragment combinations compatible with AlphaPulldown.

    Parameters:
        - proteins (list): A list of Protein objects.
        - output (str, optional): The name of the output file. Default is 'pulldown_input.txt'.
        - method (str, optional): Which protein combinations to use - 'all' for all
          v all, 'one' for all v one, 'specific' for specific combinations
          (specified using the combinations_csv argument). Default is 'all'.
        - combinations_csv (csv path, optional): The path to a CSV file containing
          specific protein combinations to generate. Each row in the CSV should
          represent one combination, with two columns for the names of the two proteins.
          Default is None. Required if method is 'specific', otherwise ignored.
        - one_protein (str, optional): Name of a protein to use for one v all
          combinations. Default is None. Required if method is 'one', otherwise
          ignored.
    
    Returns:
        - None, but creates a txt file with protein fragment combinations.
    """
    protein_pairs = get_protein_combinations(proteins, method, combinations_csv, one_protein)
    
    pulldown_lines = []
    for protein1, protein2 in protein_pairs:
        for i, (start1, end1) in enumerate(protein1.fragment_list):
            for j, (start2, end2) in enumerate(protein2.fragment_list):
                line = f"{protein1.accession_id},{start1+1}-{end1+1};{protein2.accession_id},{start2+1}-{end2+1}"
                pulldown_lines.append(line)
    with open(output_name, 'w') as file:
        file.write('\n'.join(pulldown_lines))
    print(f"File created: pulldown_combinations.txt")