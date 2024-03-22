import pandas as pd
from classes import Protein  # Importing the Protein class
from UniProt_data_fetch import fetch_uniprot_info  # Importing the function to fetch sequence data

def initialize_proteins_from_csv(csv_path):
    """
    Reads a CSV file containing protein names and accession IDs, attempts to fetch their sequences,
    and initializes a list of Protein objects. Proteins for which sequences cannot be fetched due to 
    errors or missing data are tracked separately.

    Parameters:
    - csv_path (str): The path to the CSV file.

    Returns:
    - tuple of (list of Protein, list of str): 
        - The first element is a list of initialized Protein objects for which sequences were successfully fetched.
        - The second element is a list of protein names for which sequences could not be retrieved or other errors occurred during fetching.

    This function attempts to fetch protein sequence data using provided accession IDs. If successful, a Protein object
    is created and added to the list of proteins. If the sequence cannot be fetched or an error occurs (e.g., invalid 
    UniProt code, network issues), the protein's name is added to a separate list indicating an error or missing data.
    This allows for easy identification of proteins that could not be processed as expected.
    """
    proteins = []  # List to store successfully initialized Protein objects
    proteins_with_errors = []  # List to store protein names with no UniProt data available

    # Read the CSV file
    df = pd.read_csv(csv_path, usecols=[0, 1])  # Adjust column names as necessary
    
    # Initialize an empty list to store Protein objects
    protein_list = []
    
    # Iterate through each row in the DataFrame
    for _, row in df.iterrows():
        protein_name = row[0]
        accession_id = row[1]
        try:
            # Attempt to fetch sequence using the accession ID
            fetch_result = fetch_uniprot_info(accession_id)
            sequence = fetch_result['sequence'] if 'sequence' in fetch_result else None

            if sequence:
                # Initialize a Protein object and add it to the list if sequence is found
                proteins.append(Protein(name=protein_name, accession_id=accession_id, sequence=sequence))
            else:
                # If no sequence is found, add the protein name to the error list
                proteins_with_errors.append(protein_name)
        except Exception as e:
            # If any error occurs during data fetching, also add to the error list
            proteins_with_errors.append(protein_name)

    # Return both the list of Protein objects and the list of protein names with errors
    return proteins, proteins_with_errors