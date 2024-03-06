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
    for key, value_list in dict.items():
        if value in value_list:
            return key  # Return the key if the value is found
    return None  # Return None if the value is not found

def find_domains_from_PAE(PAE):
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


# check if res1 is an item in domain_dict
# check if res2 is an item in domain_dict
# if res1 in domain_dict and res2 not in domain_dict:
    # get list of values associated with same key as res1
    # find smallest value in this list
    # make list of all values between this smallest value and res2
    # add all of these values to the items associated with the same key as res1
# if neither res1 or res2 is item in domain_dict:
    # add new key to domain_dict - key = "D+{next_domain}" - make sure this value will not change when next_domain is updated
    # associate res1, res2 and all values in between with this key (in value order from small to large - res2 will always be smaller than res1)
    # add one to next_domain



