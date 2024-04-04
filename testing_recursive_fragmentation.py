from AFprep_func.alphafold_db_domain_identification import read_afdb_json
from AFprep_func.alphafold_db_domain_identification import find_domains_from_pae
from AFprep_func.protein_init import initialize_proteins_from_csv  # Importing the function to fetch sequence data
from AFprep_func.fragment_protein import fragment_protein
from AFprep_func.long_domains import produce_and_compile_fragments

csv_path = 'AFprep_func/Asl.csv'
#csv_path = 'AFprep_func/proteins1.csv'
proteins, errors = initialize_proteins_from_csv(csv_path)
print(f"Successfully initialized proteins: {[protein.name for protein in proteins]}")
print(f"Proteins with errors or no data available: {errors}")

for protein in proteins:
    print(protein)

for protein in proteins:
    protein_PAE = read_afdb_json(protein.accession_id)
    domains = find_domains_from_pae(protein_PAE)
    for domain in domains:
        print(domain)
        protein.add_domain(domain)
    print(protein.name, " has ", len(protein.domain_list), "domains")

for protein in proteins:
    print(len(protein.sequence))
    print(protein.first_res)

for protein in proteins:
    #fragments = fragment_protein(protein)
    #for fragment in fragments:
    #    protein.add_fragment(fragment[0], fragment[1])
    produce_and_compile_fragments(protein)
    print(protein.name, " is ", len(protein.sequence), " residues long and has ", len(protein.fragment_list), "fragments")



