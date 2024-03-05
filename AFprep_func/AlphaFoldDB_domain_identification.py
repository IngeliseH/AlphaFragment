#identifying domains from AlphaFold DB
#importing required packages
import os
import requests

def read_AFDB_json(accession_id, database_version="v4"):
    """
    Fetches and returns the JSON content from the AlphaFold Database (AFDB) for a given accession ID.

    This function queries the AFDB for the predicted aligned error file associated with a specific accession ID. It handles errors where a file associated to the accession ID cannot be retrieved from the AFDB.

    Parameters:
    - accession_id (str): The accession ID for which to fetch the corresponding AlphaFold predicted aligned error JSON.
    - database_version (str, optional): The version of the database to query. Default is "v4".

    Returns:
    - dict or None: The JSON content as a dictionary if the request is successful; otherwise, None.

    Prints an error message and returns None if:
    - The accession ID is not found in the predefined dictionary.
    - The HTTP request to retrieve the file fails (e.g., file not found, network problems).

    Note:
    - This function requires the `requests` library to make HTTP requests.
    - The function assumes that alphafold ids take the form AF-[a UniProt accession]-F1. If there are multiple fragments associated with a uniprot id this will only take fragment 1
    """
    alphafold_id = f'AF-{accession_id}-F1'
    
    json_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_id}-predicted_aligned_error_{database_version}.json'
    
    try:
        response = requests.get(json_url)
        response.raise_for_status()  # Raises an HTTPError if the status is 4xx, 5xx
    except requests.exceptions.HTTPError as e:
        # Specific handling for HTTP errors (e.g., file not found on the server)
        print(f"HTTP Error: Could not retrieve file from {json_url}. Error message: {e}")
        return None
    except requests.exceptions.RequestException as e:
        # Broad exception for other issues, like network problems
        print(f"Error: A problem occurred when trying to retrieve {json_url}. Error message: {e}")
        return None
    
    # If everything went smoothly, return the content of the file
    return response.json()

#testing using Ana2
read_AFDB_json("Q9XZ31")