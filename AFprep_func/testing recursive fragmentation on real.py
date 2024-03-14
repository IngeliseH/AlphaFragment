#testing recursive fragmentation on real protein
from recursive_fragment_finding import recursive_fragmentation
from classes import Domain
from classes import Protein
from AlphaFoldDB_domain_identification import read_AFDB_json
from AlphaFoldDB_domain_identification import find_domains_from_PAE
from protein_init import initialize_proteins_from_csv  # Importing the function to fetch sequence data

# Example usage
csv_path = 'AFprep_func/proteins2.csv'
proteins, errors = initialize_proteins_from_csv(csv_path)
print(f"Successfully initialized proteins: {[protein.name for protein in proteins]}")
print(f"Proteins with errors or no data available: {errors}")

for protein in proteins:
    protein_PAE = read_AFDB_json(protein.accession_id)
    domains = find_domains_from_PAE(protein_PAE)
    for domain in domains:
        print(domain)
        protein.add_domain(domain)
    print(protein.name, " has ", len(protein.domain_list), "domains")

min_len = 100
max_len = 200
overlap = 10
max_overlap = 30
min_overlap = 0
for protein in proteins:
    if not isinstance(protein, Protein):
        raise ValueError("domain must be an instance of Domain.")
    fragments = None
    first_res = 0
    last_res = len(protein.sequence)
    while fragments == None and min_len > 0 and max_len < len(protein.sequence):
        fragments = recursive_fragmentation(protein.domain_list, first_res, first_res, last_res, min_len, max_len, overlap, max_overlap, min_overlap)
        if fragments == None:
            min_len+=-10
            print(min_len)
            max_len+=10
            print(max_len)
    if fragments is not None:
        print("Cutpoints found for protein ", protein, ":", fragments)
    else:
        print("No valid fragmentation pattern found for protein ", protein, ":")
