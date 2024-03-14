from AlphaFoldDB_domain_identification import read_AFDB_json
from AlphaFoldDB_domain_identification import find_domains_from_PAE
from protein_init import initialize_proteins_from_csv  # Importing the function to fetch sequence data
from classes import Domain

csv_path = 'proteins1.csv'
proteins, errors = initialize_proteins_from_csv(csv_path)
print(f"Successfully initialized proteins: {[protein.name for protein in proteins]}")
print(f"Proteins with errors or no data available: {errors}")

for protein in proteins:
    protein_PAE = read_AFDB_json(protein.accession_id)
    domains = find_domains_from_PAE(protein_PAE)
    for domain in domains:
        protein.add_domain(domain)
    print(protein.name, " has ", len(protein.domain_list), "domains")
