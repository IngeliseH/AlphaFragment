import requests

def fetch_uniprot_info(accession_id):
    """
    Fetches and returns protein domain information from the UniProt database for a given accession code.

    This function queries the European Bioinformatics Institute (EBI) Protein Features API for domain information associated with a specific protein accession code. It parses and returns the JSON response containing detailed information about the protein's features, such as domains, motifs, and binding sites.

    Parameters:
    - accession_id (str): The accession id of the protein for which to fetch domain information.

    Returns:
    - dict or None: The JSON response as a dictionary if the request is successful; otherwise, None.

    Prints an error message and returns None if:
    - The HTTP request fails (e.g., accession code not found, network issues).
    - An exception occurs during the request or JSON parsing.

    Note:
    - This function requires the `requests` library to make HTTP requests.
    """
    requestURL = f"https://www.ebi.ac.uk/proteins/api/features/{accession_id}"
    try:
        response = requests.get(requestURL, headers={"Accept": "application/json"})
        response.raise_for_status()  # Raises an HTTPError for bad responses (4xx, 5xx)
        return response.json()  # Parse and return JSON response
    except requests.exceptions.HTTPError as e:
        # Specific handling for HTTP errors, like 404 or 500
        print(f"HTTP Error: Could not retrieve data for accession code {accession_id}. Error message: {e}")
    except requests.exceptions.RequestException as e:
        # Broad exception for other request issues, such as network problems
        print(f"Error: A problem occurred when trying to fetch data for {accession_id}. Error message: {e}")
    return None