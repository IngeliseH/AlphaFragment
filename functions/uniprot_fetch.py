"""
This module provides functionalities for fetching and processing protein
information from the UniProt database.

Functions:
  - fetch_uniprot_info(accession_id): Fetches information for a given protein
    accession code from the UniProt database.
  - fetch_uniprot_domains(protein): Identifies and returns the domains of a
    protein based on the fetched UniProt data

Dependencies: 
  - requests: Required for making HTTP requests to the UniProt API.
  - .classes: Contains the Domain class used to represent protein
    domains.
"""

import requests
from .classes import Domain

def fetch_uniprot_info(accession_id):
    """
    Fetches and returns protein information from the UniProt database for
    a given accession code.

    Parameters:
      - accession_id (str): The accession id of the protein for which to fetch
        domain information. If 'na' or a similar placeholder is provided, the
        function will not attempt a request and will return None immediately.

    Returns:
      - dict or None: The JSON response as a dictionary if the request is
        successful; otherwise, None.

    Errors and Exceptions:
      - If the HTTP request fails due to reasons such as an invalid accession
        code or network issues, an error message is printed, and the function
        returns None.
      - If exceptions occur  during the request or while parsing the JSON
        response (eg timeouts, connection errors, or issues with the format of
        the received data), an error message is printed, and the function returns
        None.

    Note:
      - This function requires the `requests` library to make HTTP requests.
      - A timeout is set to 30 seconds for the HTTP request to prevent hanging.
    """
    # Check for non-applicable accession_id before attempting the request
    if not accession_id or accession_id.lower() == "na":
        print("No valid accession ID provided. Skipping UniProt fetch operation.")
        return None

    request_url = f"https://www.ebi.ac.uk/proteins/api/features/{accession_id}"
    try:
        response = requests.get(request_url, headers={"Accept": "application/json"}, timeout=30)
        response.raise_for_status()  # Raises an HTTPError for bad responses (4xx, 5xx)
        return response.json()  # Parse and return JSON response
    except requests.exceptions.HTTPError as e:
        # Specific handling for HTTP errors, like 404 or 500
        print(f"HTTP Error: Could not retrieve data for accession code "
              f"{accession_id}. Error message: {e}")
    except requests.exceptions.Timeout:
        print(f"Timeout Error: The request timed out while attempting to fetch "
              f"data for {accession_id}.")
    except requests.exceptions.RequestException as e:
        # Broad exception for other request issues, such as network problems
        print(f"Error: A problem occurred when trying to fetch data for "
              f"{accession_id}. Error message: {e}")
    return None

def find_uniprot_domains(protein):
    """
    Extracts domain information for a given protein from UniProt data.

    Parameters:
      - protein (Protein): The protein object for which domain information is to
        be fetched. The object should have a valid 'accession_id' attribute that
        corresponds to its UniProt accession code.

    Returns:
      - list of Domain objects: A list containing Domain objects for each domain
        found in the UniProt data. Domain positions are given using Pythonic 0-based
        indexing, so will be 1 less than the UniProt positions.

    Errors and Exceptions:
      - If no data can be retrieved from UniProt for by the `fetch_uniprot_info`
        function returns None, this function prints a message and returns None.
      - If the start or end positions of a domain cannot be converted to integers,
        a ValueError will be raised.

    Notes:
      - The function intentionally excludes domains of types 'CHAIN', 'MUTAGEN', and
        'CONFLICT'. This is because:
        - CHAIN types cover the entire sequence and do not represent discrete structural
          domains, making them unsuitable for analyses that require domain fragmentation.
        - MUTAGEN and CONFLICT types are related to sequence variations and conflicts,
          respectively, and are not considered structurally relevant domains for the
          purpose of protein analysis and domain identification.
      - In specific cases where it is necessary to include these excluded domain types,
        or to specify domains not found or incorrectly annotated in UniProt, users
        can manually specify domains in the input CSV file used with other functions of
        this module.
    """
    data = fetch_uniprot_info(protein.accession_id)
    uniprot_domains = []

    if data is None:
        return None

    for feature in data['features']:
        if feature['type'] not in ['CHAIN', 'MUTAGEN', 'CONFLICT']:
            description = feature.get('description', 'No description')
            try:
                start = int(feature['begin'])
                end = int(feature['end'])
            except ValueError:
                print(f"Invalid start/end in feature for protein {protein.name}: {feature}")
                continue

            uniprot_domains.append(Domain(description, start-1, end-1, "UniProt"))
    if uniprot_domains:
        print(f"{len(uniprot_domains)} domains found in UniProt for "
              f"{protein.name}: {uniprot_domains}")
    else:
        print(f"No domains found in UniProt for protein {protein.name}.")
    return uniprot_domains
