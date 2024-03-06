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


# Function to find the key for a given value
def find_key_for_value(dict, value):
    """
    Searches through a dictionary to find the key associated with a given value.

    This function iterates over each key-value pair in the dictionary. It checks whether the given value is present
    in the list of values associated with each key. If the value is found within a value list, the function returns
    the corresponding key.

    Parameters:
    - dict (dict): The dictionary to search through. It is expected that the dictionary's values are lists.
    - value: The value to search for within the lists associated with each key in the dictionary.

    Returns:
    - The key associated with the given value if found within any of the value lists. If the value is not found,
      the function returns None.

    Note:
    - This function assumes that each value in the dictionary is a list of values. If the given value is found
      within any of these lists, the corresponding key is returned.
    - If multiple keys contain the value in their respective lists, only the first key encountered during the
      iteration is returned. The order of iteration is determined by the insertion order of the keys.
    - The function returns None if the value is not found in any of the lists within the dictionary.
    """
    for key, value_list in dict.items():
        if value in value_list:
            return key  # Return the key if the value is found
    return None  # Return None if the value is not found


def find_domains_from_PAE(PAE):
    """
    Analyzes Predicted Aligned Error (PAE) data to group residues into domains based on their PAE values and distances.

    This function iterates through pairs of residues in the provided PAE matrix. It determines whether each pair belongs
    to the same domain based on their PAE values relative to specified cutoffs for "close" and "further" distances. If a pair
    is deemed to be in the same domain, they are grouped together in the `domain_dict` dictionary, with each domain
    assigned a unique identifier.

    Parameters:
    - PAE (list of lists): A 2D matrix representing the Predicted Aligned Error values between residue pairs. PAE[i][j]
      gives the PAE value between residue i and residue j.

    Returns:
    - domain_dict (dict): A dictionary where each key is a domain identifier (e.g., "D1", "D2", ...) and each value is a
      list of residue indices belonging to that domain. The indices in each list are unique and sorted in ascending order.

    The function employs a specific logic to decide domain membership:
    - Two residues are considered to be in the same domain if the distance between them is greater than `res_dist_cutoff` and the
      lesser of the two PAE values between them (ie min(PAE[res1, res2] and PAE[res2, res1])) is less than `further_PAE_val`, or if
      the distance between them is less than or equal to `res_dist_cutoff` and the lesser of the two PAE values between them is below
      `close_PAE_val`. A higher PAE threshold is set for closer residues, as these will always have a higher background level of
      confidence about their relative positions.
    - PAE between residues with less than 4 residues between them will not be examined, due to always having high confidence in
      this range, even in unstructured regions
    - If either residue in the pair is already associated with a domain, the pair is added to the existing domain. If both are
      new, a new domain is created. If only one residue is new, it is added to the domain of the other residue.
    - The inner loop breaks early once a residue pair is processed and determined to be in the same domain, moving to the next residue.

    Note:
    - `res_dist_cutoff`, `close_PAE_val`, and `further_PAE_val` are predefined thresholds used to evaluate the PAE data.
    - The function assumes that the PAE matrix is symmetric and that the PAE value between a residue and itself is not considered.
    """
    domain_dict = {}
    next_domain = 1
    res_dist_cutoff = 7
    close_PAE_val = 2
    further_PAE_val = 6

    for res1 in range(len(PAE)-1, -1, -1):
        for res2 in range(0, res1 - 3):
            # Calculate the distance between resiudes being evaluated
            res_difference = abs(res1 - res2)
            # Find the PAE between the residues, looking at both directions
            relative_PAE = min(PAE[res1, res2], PAE[res2, res1])

            # Evaluate whether residues are part of the same domain given PAE value and their distance
            is_same_domain = (res_difference <= res_dist_cutoff and relative_PAE < close_PAE_val) or \
                             (res_difference > res_dist_cutoff and relative_PAE < further_PAE_val)

            if is_same_domain:
                key_res1 = find_key_for_value(domain_dict, res1)
                key_res2 = find_key_for_value(domain_dict, res2)

                if key_res1 and not key_res2:
                    # Add res2 and all values in between to the domain of res1
                    domain_dict[key_res1].extend(range(min(domain_dict[key_res1]), res2 + 1))
                    domain_dict[key_res1] = list(set(domain_dict[key_res1]))  # Remove duplicates
                elif not key_res1 and not key_res2:
                    # Create new domain and associate res1, res2, and all values in between
                    domain_dict[f"D{next_domain}"] = list(range(res2, res1 + 1))
                    next_domain += 1

                break  # Break the inner loop once is_same_domain condition is met and processed
    return domain_dict    

#testing using Ana2
Ana2_domains = find_domains_from_PAE(Ana2_PAE)
print(Ana2_domains)





