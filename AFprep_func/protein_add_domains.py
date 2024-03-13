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
    domain_dict = find_domains_from_PAE(protein_PAE)
    for domain_name, residues in domain_dict.items():
        start = min(residues)
        end = max(residues)
        domain_type = 'AF'
        domain = Domain(start=start, end=end, domain_type=domain_type, residues=residues)
        protein.add_domain(domain)
    print(protein.name, " has ", len(protein.domain_list), "domains")
