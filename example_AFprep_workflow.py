from AFprep_func.alphafold_db_domain_identification import read_afdb_json
from AFprep_func.alphafold_db_domain_identification import find_domains_from_pae
from AFprep_func.protein_init import initialize_proteins_from_csv
from AFprep_func.fragment_protein import fragment_protein
from AFprep_func.plot_fragments import plot_fragmentation_output

#initialising proteins from csv
CSV_PATH = 'AFprep_func/Asl.csv'
#CSV_PATH = 'AFprep_func/proteins1.csv'
proteins, errors = initialize_proteins_from_csv(CSV_PATH)
print(f"Successfully initialized proteins: {[protein.name for protein in proteins]}")
print(f"Proteins with errors or no data available: {errors}")

#finding domains from AlphaFold
for protein in proteins:
    protein_PAE = read_afdb_json(protein.accession_id)
    domains = find_domains_from_pae(protein_PAE)
    for domain in domains:
        print(domain)
        protein.add_domain(domain)
    print(protein.name, " has ", len(protein.domain_list), "domains")

#fragmenting proteins
for protein in proteins:
    fragments = fragment_protein(protein)
    print("Cutpoints found for protein ", protein.name, ":", fragments)
    for start, end in fragments:
        protein.add_fragment(start, end)
    plot_fragmentation_output(protein, protein.fragment_list, save_location = None)
    print(protein.name, " is ", len(protein.sequence), " residues long and has ",
          len(protein.fragment_list), "fragments")
