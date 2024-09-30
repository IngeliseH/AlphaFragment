"""
This module contains functions for creating protein fragment pairs and files using
these, for use as input to AlphaFold or AlphaPulldown.

Functions:
    - get_protein_combinations: Generates a list of protein pairs based on the method selected.
    - output_fastas: Creates fasta files for each combination of protein fragments.
    - output_pulldown: Creates a txt file with protein fragment combinations compatible
      with AlphaPulldown.

Dependencies:
    - csv: For reading CSV files.
    - os: For creating directories and files.
    - combinations: For generating combinations of proteins.
    - Protein: A class representing a protein with sequence and domain information.
"""
import csv
import os
from itertools import combinations
from .classes import Protein

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
    
    Note:
        - For 'all' and 'one' methods, self-combinations are included.
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

    # For 'all' method, compile all combinations of proteins
    if method == 'all':
        protein_pairs = list(combinations(proteins, 2)) + [(p, p) for p in proteins]

    # For 'one' method, find target protein + compile all combinations with it
    elif method == 'one':
        target_protein = next((p for p in proteins if p.name == one_protein), None)
        if not target_protein:
            raise ValueError(f"Protein named {one_protein} not found among the provided proteins.")
        protein_pairs = [(target_protein, p) for p in proteins]

    # For 'specific' method, read the desired combinations from the CSV
    elif method == 'specific':
        protein_dict = {protein.name: protein for protein in proteins}
        protein_pairs = []
        with open(combinations_csv, mode='r', encoding='utf-8', newline='') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                if row[0] and row[1]:
                # Check if both proteins specified in the row exist in the protein list
                    if row[0] in protein_dict and row[1] in protein_dict:
                        protein_pairs.append((protein_dict[row[0]], protein_dict[row[1]]))
                    else:
                        # Print error message if one or both proteins are not found
                        missing_proteins = [name for name in row[:2] if name not in protein_dict]
                        print(f"Error: Combination {row[0]}-{row[1]} not possible. "
                              f"Missing proteins: {', '.join(missing_proteins)}.")

    return sorted(protein_pairs, key=lambda pair: (pair[0].name, pair[1].name))

def output_fastas(proteins, save_location=None, method='all', combinations_csv=None, one_protein=None):
    """
    Creates a folder for each protein pair, containing fasta files for each combination
    of fragments for that pair.
    
    Parameters:
        - proteins (list): A list of Protein objects.
        - save_location (str, optional): The path to the directory where the folders
          and fasta files will be saved. Default is None, which saves in the current
          working directory.
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
        if save_location:
            folder_path = os.path.join(save_location, f"{protein1.name}_{protein2.name}")
        else:
            folder_path = f"{protein1.name}_{protein2.name}"
        os.makedirs(folder_path, exist_ok=True)

        fragment_pairs=[]
        for i, (start1, end1) in enumerate(protein1.fragment_list):
            for j, (start2, end2) in enumerate(protein2.fragment_list):
                #don't create duplicate files (fragment1, fragment2) == (fragment2, fragment1)
                if (start1, end1, start2, end2) in fragment_pairs or (
                    start2, end2, start1, end1) in fragment_pairs:
                    continue

                #if not duplicate, create file
                fragment_pairs.append((start1, end1, start2, end2))
                filename = os.path.join(folder_path,
                                        f"{protein1.name}_F{i+1}+{protein2.name}_F{j+1}.fasta")
                content = (f">{protein1.name}_F{i+1}+{protein2.name}_F{j+1}\n"
                           f"{protein1.sequence[start1:end1]}:{protein2.sequence[start2:end2]}")
                with open(filename, 'w', encoding='utf-8') as file:
                    file.write(content)
                print(f"File created: {filename}")

def output_pulldown(proteins, output_name='pulldown_input.txt', method='all',
                    combinations_csv=None, one_protein=None):
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
        fragment_pairs=[]
        for _, (start1, end1) in enumerate(protein1.fragment_list):
            for _, (start2, end2) in enumerate(protein2.fragment_list):
                # don't create duplicate files (fragment1, fragment2) == (fragment2, fragment1)
                if (start1, end1, start2, end2) in fragment_pairs or (
                    start2, end2, start1, end1) in fragment_pairs:
                    continue
                line = (f"{protein1.accession_id},{start1+1}-{end1};"
                        f"{protein2.accession_id},{start2+1}-{end2}")
                pulldown_lines.append(line)
                fragment_pairs.append((start1, end1, start2, end2))

    with open(output_name, 'w',  encoding='utf-8') as file:
        file.write('\n'.join(pulldown_lines))
    print(f"File created: {output_name}")
