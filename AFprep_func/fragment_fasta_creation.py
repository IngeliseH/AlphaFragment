"""
This module contains functions for creating fasta files of protein fragment pairs,
for use as input to AlphaFold

Functions:
    - read_protein_combinations: Reads protein combinations from a CSV file.
    - generate_protein_combination_files: Generates folders and .fasta files for
      combinations of protein fragments.

Dependencies:
    - csv: For reading CSV files.
    - os: For creating directories and files.
    - combinations: For generating combinations of proteins.
"""
import csv
import os
from itertools import combinations
from AFprep_func.classes import Protein

def read_protein_combinations(csv_path):
    """
    Reads protein combinations from a CSV file.

    Parameters:
        - csv_path (str): The path to the CSV file containing protein combinations.
          Each row in the CSV should represent one combination, with two columns
          for the names of the two proteins.
    
    Returns:
        - protein_combinations (list): A list of tuples, where each tuple contains
          the names of two proteins that should be combined.
    """
    protein_combinations = []
    with open(csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if row[0] and row[1]:
                protein_combinations.append((row[0], row[1]))
    return protein_combinations

def generate_protein_combination_files(proteins, method=all, combinations_csv=None, one_protein=None):
    """
    Generates a folder for each protein pair, containing a .fasta file for each
    combination of protein fragments within that protein pair.

    Parameters:
        - proteins (list): A list of Protein objects.
        - method (str, optional): which protein combinations to use - 'all' for all
          v all, 'one' for all v one, 'specific' for specific combinations
          (specified using the specific_combinations argument). Default is 'all'.
        - combinations_csv (csv path, optional): The path to a CSV file containing
          specific protein combinations to generate. Each row in the CSV should
          represent one combination, with two columns for the names of the two proteins.
          Default is None. Required if method is 'specific', otherwise ignored.
        - one_protein (str, optional): Name of a protein to use for one v all
          combinations. Default is None. Required if method is 'one', otherwise
          ignored.
    
    Returns:
        - None, but creates folders and .fasta files for each protein pair.
    
    Raises:
        - ValueError: If the input proteins are not Protein objects, if specific
          combinations are required but not provided, or if a target protein is
          required but not provided.
    """
    # Check that proteins is a list of Protein objects
    if not all(isinstance(protein, Protein) for protein in proteins):
        raise ValueError("Input proteins must be Protein objects.")
    # Check that the method is valid
    if method not in ['all', 'one', 'specific']:
        raise ValueError("Method must be 'all', 'one', or 'specific'.")
    # Check that specific combinations are provided if required
    if method == 'specific' and combinations_csv is None:
        raise ValueError("Specific combinations selected but no combinations provided.")
    # Check that a target protein is provided if required
    if method == 'one' and not one_protein:
        raise ValueError("Method 'one' selected but no target protein specified.")

    # Read in protein combinations if using specific combinations
    if method == 'specific':
        specific_combinations = read_protein_combinations(combinations_csv)

    # If using one v all, find the target protein in proteins list
    target_protein = None
    if method == 'one':
        for protein in proteins:
            if protein.name == one_protein:
                target_protein = protein
                break
        if not target_protein:
            raise ValueError(f"Protein named {one_protein} not found among the provided proteins.")

    # Pair every protein with every other, including itself
    protein_pairs = list(combinations(proteins, 2)) + [(protein, protein) for protein in proteins]

    for pair in protein_pairs:
        protein1, protein2 = pair
        # Filter combinations according to the selected method
        if method == 'specific' and ((protein1.name, protein2.name) not in specific_combinations and (protein2.name, protein1.name) not in specific_combinations):
            continue
        if method == 'one' and (protein1 != target_protein and protein2 != target_protein):
            continue

        folder_name = f"{protein1.name}_{protein2.name}"
        os.makedirs(folder_name, exist_ok=True)

        for i, (start1, end1) in enumerate(protein1.fragment_list, start=1):
            fragment1 = protein1.sequence[start1:end1]
            for j, (start2, end2) in enumerate(protein2.fragment_list, start=1):
                fragment2 = protein2.sequence[start2:end2]
                filename = os.path.join(folder_name, f"{protein1.name}_F{i}+{protein2.name}_F{j}.fasta")
                content = f">{protein1.name}_F{i}+{protein2.name}_F{j}\n{fragment1}:{fragment2}"
                with open(filename, 'w') as file:
                    file.write(content)
                print(f"File created: {filename}")