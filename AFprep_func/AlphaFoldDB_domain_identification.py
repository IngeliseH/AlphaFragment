#identifying domains from AlphaFold DB
#importing required packages
import os
import requests

def read_AFDB_json(accession_id, database_version="v4"):
    """
    Fetches and returns the Predicted Aligned Error (PAE) data from the AlphaFold Database (AFDB) for a given UniProt accession ID.

    This function queries the AFDB for the PAE associated with a specific accession ID and the chosen database version. It directly extracts and returns the 'predicted_aligned_error' data from the JSON response.

    Parameters:
    - accession_id (str): The accession ID for which to fetch the corresponding AlphaFold PAE data.
    - database_version (str, optional): The version of the database to query. Default is "v4".

    Returns:
    - list of lists or None: The predicted_aligned_error data as a list of lists if the request is successful and the data is present; otherwise, None.

    Prints an error message and returns None if:
    - The HTTP request to retrieve the file fails (e.g., file not found, network problems).
    - The 'predicted_aligned_error' data cannot be found within the response, indicating either an issue with the response data structure or the absence of PAE data for the provided accession ID.


    Note:
    - This function requires the `requests` library to make HTTP requests.
    - The function assumes that alphafold ids take the form AF-[a UniProt accession]-F1. If there are multiple fragments associated with a uniprot id this will only take fragment 1
    """
    alphafold_id = f'AF-{accession_id}-F1'
    json_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_id}-predicted_aligned_error_{database_version}.json'
    
    try:
        response = requests.get(json_url)
        response.raise_for_status()  # Raises an HTTPError if the status is 4xx, 5xx
        data = response.json()
        # Extract just the 'predicted_aligned_error' data
        predicted_aligned_error = data[0]["predicted_aligned_error"] if data and "predicted_aligned_error" in data[0] else None
    except requests.exceptions.HTTPError as e:
        # Specific handling for HTTP errors (e.g., file not found on the server)
        print(f"HTTP Error: Could not retrieve file from {json_url}. Error message: {e}")
        return None
    except requests.exceptions.RequestException as e:
        # Broad exception for other issues, like network problems
        print(f"Error: A problem occurred when trying to retrieve {json_url}. Error message: {e}")
        return None
    
    # If everything went smoothly, return the content of the file
    return predicted_aligned_error

#testing using Ana2
Ana2_PAE = read_AFDB_json("Q9XZ31")
print(Ana2_PAE)

