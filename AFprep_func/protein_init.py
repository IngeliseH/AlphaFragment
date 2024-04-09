"""
This module provides functionality to turn a CSV file listing protein names
and associated accession ids into a list of Protein objects, using sequence
information from UniProt. Proteins for which sequence data cannot be fetched are
tracked and reported separately.

Functions:
- initialize_proteins_from_csv: Reads a CSV file and initializes Protein objects
  with sequence data from UniProt.

Dependencies:
- pandas: Used for reading and processing the CSV file.
- AFprep_func.classes.Protein: The Protein class used to represent and manipulate
  protein data.
- AFprep_func.uniprot_data_fetch.fetch_uniprot_info: A function to fetch protein
  data from UniProt.

This module is intended for bioinformatics applications, specifically within the
AlphaPrep workflow
"""

import pandas as pd
from requests.exceptions import RequestException
# Import the Protein class
from AFprep_func.classes import Protein
# Import function to fetch sequence data
from AFprep_func.uniprot_data_fetch import fetch_uniprot_info

def initialize_proteins_from_csv(csv_path):
    """
    Reads a CSV file containing protein names and accession IDs, attempts to
    fetch their sequences, and initializes a list of Protein objects. Proteins
    for which sequences cannot be fetched due to errors or missing data are
    tracked separately.

    Parameters:
    - csv_path (str): The path to the CSV file.

    Returns:
    - tuple of (list of Protein, list of str): 
        - The first element is a list of initialized Protein objects for which
          sequences were successfully fetched.
        - The second element is a list of protein names for which sequences
          could not be retrieved or other errors occurred during fetching.

    This function attempts to fetch protein sequence data using provided
    accession IDs. If successful, a Protein object is created and added to the
    list of proteins. If the sequence cannot be fetched or an error occurs (e.g.,
    invalid UniProt code, network issues), the protein's name is added to a
    separate list indicating an error or missing data. This allows for easy
    identification of proteins that could not be processed as expected.
    """
    proteins = []  # List to store successfully initialized Protein objects
    proteins_with_errors = []  # List to store protein names with no UniProt data available

    # Read the CSV file
    df = pd.read_csv(csv_path, usecols=[0, 1], names=['name', 'accession_id'])

    # Iterate through each row in the DataFrame
    for _, row in df.iterrows():
        protein_name = row['name']
        print("protein name = ", protein_name)
        accession_id = row['accession_id']
        print("accession id = ", accession_id)
        try:
            # Attempt to fetch sequence using the accession ID
            fetch_result = fetch_uniprot_info(accession_id)
            if fetch_result is not None and 'sequence' in fetch_result:
                sequence = fetch_result['sequence']
                # Initialize a Protein object and add it to the list if sequence is found
                proteins.append(Protein(name=protein_name,
                                        accession_id=accession_id,
                                        sequence=sequence))
            else:
                # If no result is found, add the protein name to the error list
                proteins_with_errors.append(protein_name)
        except RequestException as e:
            print(f"Error fetching data for {protein_name}: {str(e)}")
            # If any error occurs during data fetching, also add to the error list
            proteins_with_errors.append(protein_name)

    # Return both the list of Protein objects and the list of protein names with errors
    return proteins, proteins_with_errors
