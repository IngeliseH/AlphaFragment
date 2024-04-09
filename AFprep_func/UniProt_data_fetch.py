"""
This module provides functionalities for fetching domain information associated
with protein accession codes from the UniProt database.

Dependencies:
  - requests: Required for making HTTP requests to the UniProt API.

Functions:
  - fetch_uniprot_info(accession_id): Fetches and returns information for a
    given protein accession code from the UniProt database.
"""

import requests

def fetch_uniprot_info(accession_id):
    """
    Fetches and returns protein information from the UniProt database for
    a given accession code.

    Parameters:
      - accession_id (str): The accession id of the protein for which to fetch
        domain information.

    Returns:
      - dict or None: The JSON response as a dictionary if the request is
        successful; otherwise, None.

    Prints an error message and returns None if:
      - The HTTP request fails (e.g., accession code not found, network issues).
      - An exception occurs during the request or JSON parsing.

    Note:
      - This function requires the `requests` library to make HTTP requests.
      - A timeout is set to 30 seconds for the HTTP request to prevent hanging.
    """
    request_url = f"https://www.ebi.ac.uk/proteins/api/features/{accession_id}"
    try:
        response = requests.get(request_url, headers={"Accept": "application/json"}, timeout=30)
        response.raise_for_status()  # Raises an HTTPError for bad responses (4xx, 5xx)
        return response.json()  # Parse and return JSON response
    except requests.exceptions.HTTPError as e:
        # Specific handling for HTTP errors, like 404 or 500
        print(f"HTTP Error: Could not retrieve data for accession code {accession_id}. Error message: {e}")
    except requests.exceptions.Timeout:
        print(f"Timeout Error: The request timed out while attempting to fetch data for {accession_id}.")
    except requests.exceptions.RequestException as e:
        # Broad exception for other request issues, such as network problems
        print(f"Error: A problem occurred when trying to fetch data for {accession_id}. Error message: {e}")
    return None
