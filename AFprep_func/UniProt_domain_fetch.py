def fetch_domain_info(accession_code):
    requestURL = f"https://www.ebi.ac.uk/proteins/api/features/{accession_code}"
    response = requests.get(requestURL, headers={"Accept": "application/json"})
    
    if response.ok:
        return response.json()  # Parse JSON response
    else:
        response.raise_for_status()