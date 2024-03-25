#identifying domains from AlphaFold DB
#importing required packages
import requests
from AFprep_func.classes import Domain

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

def find_domains_from_PAE(PAE):
    """
    Analyzes Predicted Aligned Error (PAE) data to group residues into domains. This function iterates through residue
    pairs, determining their domain membership based on PAE values and residue distances. Domains are represented as 
    Domain objects with unique identifiers, start, and end residues.

    Parameters:
    - PAE (list of lists): A 2D matrix of PAE values between residue pairs, where PAE[i][j] is the PAE between residues i and j.

    Returns:
    - A list of Domain objects, each representing a domain with a unique identifier and the range of residues it encompasses.

    The logic used to determine domain membership is:
    
    - Two residues are considered to be in the same domain if the distance between them is greater than `res_dist_cutoff` and the
      lesser of the two PAE values between them (ie min(PAE[res1, res2] and PAE[res2, res1])) is less than `further_PAE_val`, or if
      the distance between them is less than or equal to `res_dist_cutoff` and the lesser of the two PAE values between them is below
      `close_PAE_val`. A higher PAE threshold is set for closer residues, as these will always have a higher background level of
      confidence about their relative positions.
    - The function does not evaluate PAE for residue pairs less than 4 residues apart due to inherently high confidence in their relative positions.
    - Domains are updated or created based on the membership of the residues being evaluated. If one residue is already in a domain and the other is
      not, the latter is added to the former's domain. New domains are created for pairs where neither residue is currently in a domain.
    - The inner loop breaks early once a residue pair is processed and determined to be in the same domain, moving to the next residue (as all
      residues between these are assumed to be in the same domain). No gaps are possible within domains

    Note:
    - `res_dist_cutoff`, `close_PAE_val`, and `further_PAE_val` are defined thresholds for evaluating PAE data.
    - Assumes the PAE matrix is symmetric and PAE[i][i] (self-comparison) is not considered.
    """
    domains = []
    next_domain_num = 1
    res_dist_cutoff = 10
    close_PAE_val = 4
    further_PAE_val = 11

    for res1 in range(0, len(PAE)):  # Iterate through residues from start to end
        for res2 in range(len(PAE)-1, res1 + 4, -1):  # Iterate through potential domain-mate residues, skipping nearby ones
            
            # Calculate the distance between residues being evaluated
            res_difference = abs(res2 - res1)
            # Find the PAE between the residues, looking at both directions
            relative_PAE = min(PAE[res1][res2], PAE[res2][res1])

            # Determine if residues are in the same domain based on PAE and distance
            is_same_domain = ((res_difference <= res_dist_cutoff and relative_PAE < close_PAE_val) or
                              (res_difference > res_dist_cutoff and relative_PAE < further_PAE_val))

            if is_same_domain:
                domain_res1 = find_domain_by_res(domains, res1)
                domain_res2 = find_domain_by_res(domains, res2)

                if domain_res1 and not domain_res2:
                    # Extend domain_res1 to include res2 (and so all residues in between) if outside current range
                    domain_res1.end = max(domain_res1.end, res2)
                elif not domain_res1 and not domain_res2:
                    # Create a new domain starting at res1 and ending at res2
                    domains.append(Domain(f"D{next_domain_num}", res1, res2, 'AF'))
                    next_domain_num += 1

                break  # Move to the next residue after processing a same domain pair
            
    return domains


