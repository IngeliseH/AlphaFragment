"""
This module interacts with the AlphaFold Database (AFDB) to obtain and analyze
Predicted Aligned Error (PAE) data for identifying protein domains. Domains are 
determined from PAE values, grouping residues with similar PAE thresholds.

Functions:
  - read_afdb_json(accession_id, database_version="v4"): Retrieves PAE data.
  - find_domain_by_res(domains, res): Helper function to identify the domain a
    residue is in.
  - find_domains_from_pae(pae): Groups residues into domains based on PAE.

Dependencies:
  - requests: Used for making HTTP requests to the AlphaFold Database.
  - AFprep_func.classes.Domain: The Domain class used to represent protein domains.
"""

#importing required packages
import requests
from AFprep_func.classes import Domain

def read_afdb_json(accession_id, database_version="v4"):
    """
    Fetches and returns the Predicted Aligned Error (PAE) data from the
    AlphaFold Database (AFDB) for a given UniProt accession ID.

    Parameters:
      - accession_id (str): The accession ID for which to fetch the corresponding
        AlphaFold PAE data. If 'na' or a similar placeholder is provided, the
        function will not attempt a request and will return None immediately.
      - database_version (str, optional): The version of the database to query.
        Default is "v4".

    Returns:
      - list of lists or None: The predicted_aligned_error data as a list of
        lists if the request is successful and the data is present; otherwise,
        None.

    Prints an error message and returns None if:
      - The HTTP request to retrieve the file fails (e.g., file not found, network
        problems).
      - The 'predicted_aligned_error' data cannot be found within the response,
        indicating either an issue with the response data structure or the absence
        of PAE data for the provided accession ID.

    Note:
      - This function requires the `requests` library to make HTTP requests.
      - The function assumes that alphafold ids take the form
        'AF-[a UniProt accession]-F1.' If there are multiple fragments associated
        with a uniprot id this will only take fragment 1
    """
    # Check for non-applicable accession_id before attempting the request
    if accession_id.lower() == "na":
        print("No valid accession ID provided. Skipping fetch operation.")
        return None
    
    alphafold_id = f'AF-{accession_id}-F1'
    json_url = (
        f'https://alphafold.ebi.ac.uk/files/{alphafold_id}'
        f'-predicted_aligned_error_{database_version}.json'
    )
    try:
        response = requests.get(json_url, timeout=30) #timeout in 30 seconds
        # Raise an HTTPError if the status is 4xx, 5xx
        response.raise_for_status()
        data = response.json()
        # Extract just the 'predicted_aligned_error' data
        predicted_aligned_error = (
            data[0]["predicted_aligned_error"] if data and "predicted_aligned_error" in data[0]
            else None
        )
    except requests.exceptions.HTTPError as e:
        # Specific handling for HTTP errors (e.g., file not found on the server)
        error_message = f"HTTP Error: Could not retrieve file from {json_url}."
        print(f"{error_message} Error message: {e}")
        return None
    except requests.exceptions.RequestException as e:
        # Broad exception for other issues, like network problems
        error_message = f"Error: A problem occurred when trying to retrieve {json_url}."
        print(f"{error_message} Error message: {e}")
        return None

    # If everything went smoothly, return the content of the file
    return predicted_aligned_error


# Function to find the domain for a given residue
def find_domain_by_res(domains, res):
    """
    Helper function to find the domain that contains a given residue.

    Parameters:
      - domains (list of Domain): The list of current domain objects.
      - res (int): The residue index to find within the domains.

    Returns:
      - Domain object if the residue is found within a domain, None otherwise.
    """
    for domain in domains:
        if domain.start <= res <= domain.end:
            return domain
    return None

def find_domains_from_pae(pae,  method='cautious', custom_params=None):
    """
    Analyzes Predicted Aligned Error (PAE) data to group residues into domains.
    This function iterates through residue pairs, determining their domain
    membership based on PAE values and residue distances. Domains are
    represented as Domain objects with unique identifiers, start, and end
    residues.

    Parameters:
      - pae (list of lists): A 2D matrix of PAE values between residue pairs,
        where pae[i][j] is the PAE between residues i and j.
      - method (str, optional): Strategy for grouping residues into domains.
        Options are:
          - 'cautious' - Groups residues into domains with moderate PAE
            thresholds, aiming to balance sensitivity and specificity.
          - 'definite' - Only groups residues into domains if there is very high
            confidence in their relative positions. Likely to produce smaller domains.
          - 'custom' - Allows specification of custom PAE thresholds via the
            custom_params argument.
      - custom_params (dict, optional): Required if method is 'custom', ignored
        otherwise. Must include:
          - 'res_dist_cutoff' (int) - The residue distance threshold to
            differentiate between close and further residue evaluations.
            Meaningless if close_pae_val == further_pae_val, or if set below 4
            as residues closer than this are ignored anyway due to high base
            confidence in their relative positions. Set to 10 for cautious
            grouping method, and irrelavant for definite grouping method as close
            and further pae_vals are equal.
          - 'close_pae_val' (int) - The PAE threshold which residue pairs must
            fall below to be considered within the same domain, if the distance
            between them is between 4 and res_dist_cutoff. Set to 4 for cautious
            grouping method and 2 for definite grouping method.
          - 'further_pae_val' (int) - The PAE threshold which residue pairs must
            fall below to be considered within the same domain, if the distance
            between them is greater than res_dist_cutoff. Set to 11 for close
            grouping method and 2 for definite grouping method.

    Returns:
      - A list of Domain objects, each representing a domain with a unique
        identifier and the range of residues it encompasses.

    Function logic:
      - Residues are grouped into the same domain based on a comparison of their
        Predicted Aligned Error (pae) values against thresholds determined by
        their proximity. If two residues are sufficiently close (based on a
        predefined or custom distance threshold), their PAE value must fall
        below a stricter, lower threshold to confirm high confidence in their 
        proximity. For further apart residues, a higher PAE threshold can be
        used to account for less background confidence in their relative positions.
      - If a residue pair is decided to be in the same domain, the function checks
        if either residue is already part of an existing domain, and either adds
        the other residue to that domain or creates a new domain that includes both
        residues. All residues in between the two being assessed will also be
        automatically included in the domain - no gaps are allowed within domains.
      - Very close residues will always have high confidence in their relative
        positions, so are not evaluated.

    Raises:
      - ValueError: If the method is 'custom' and custom_params is not provided
        or missing necessary keys, if an invalid method is specified, or if the
        pae data is in the wrong format (not a square matrix of numbers).
    """
    # Check if 'pae' is a non-empty matrix, each row is a list, and contains only numeric entries
    if not pae or not isinstance(pae, list) or any(not isinstance(row, list) or any(not isinstance(item, (int, float)) for item in row) for row in pae):
        raise ValueError("Input 'pae' must be a non-empty matrix of numbers.")

    # Check if 'pae' is a square matrix
    if any(len(row) != len(pae) for row in pae):
        raise ValueError("Input 'pae' must be a square matrix.")

    # Default parameters for cautious and definite methods
    parameters = {
        'cautious': {'res_dist_cutoff': 10, 'close_pae_val': 4, 'further_pae_val': 11},
        'definite': {'res_dist_cutoff': 0, 'close_pae_val': 2, 'further_pae_val': 2}
    }

    if method in ['cautious', 'definite']:
        res_dist_cutoff = parameters[method]['res_dist_cutoff']
        close_pae_val = parameters[method]['close_pae_val']
        further_pae_val = parameters[method]['further_pae_val']
    elif method == 'custom':
        if not custom_params or not all(key in custom_params for key in ['res_dist_cutoff', 'close_pae_val', 'further_pae_val']):
            raise ValueError("For custom method, 'custom_params' must be a dictionary with keys 'res_dist_cutoff', 'close_pae_val', 'further_pae_val'.")
        res_dist_cutoff = custom_params['res_dist_cutoff']
        close_pae_val = custom_params['close_pae_val']
        further_pae_val = custom_params['further_pae_val']
    else:
        raise ValueError("Invalid method. Choose 'cautious', 'definite', or 'custom'.")

    domains = []
    next_domain_num = 1

    # Iterate through residues from start to end
    for res1 in range(0, len(pae)):
        # Iterate through potential domain-mate residues, skipping nearby ones
        # Start at end of protein and work in reverse
        for res2 in range(len(pae)-1, res1 + 4, -1):

            # Calculate the distance between residues being evaluated
            res_difference = abs(res2 - res1)
            # Find the PAE between the residues, looking at both directions
            relative_pae = min(pae[res1][res2], pae[res2][res1])

            # Determine if residues are in the same domain based on PAE and distance
            is_same_domain = ((res_difference <= res_dist_cutoff and
                               relative_pae < close_pae_val) or
                              (res_difference > res_dist_cutoff and
                               relative_pae < further_pae_val))

            if is_same_domain:
                domain_res1 = find_domain_by_res(domains, res1)
                domain_res2 = find_domain_by_res(domains, res2)

                if domain_res1 and not domain_res2:
                    # Extend domain_res1 to include all res up to and including res2
                    domain_res1.end = max(domain_res1.end, res2)
                elif not domain_res1 and not domain_res2:
                    # Create a new domain starting at res1 and ending at res2
                    domains.append(Domain(f"AF_D{next_domain_num}", res1, res2, 'AF'))
                    next_domain_num += 1

                # Move to the next residue after processing a same domain pair
                break

    return domains
