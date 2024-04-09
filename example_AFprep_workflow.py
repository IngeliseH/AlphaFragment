import pandas as pd
from ast import literal_eval
from AFprep_func.protein_init import initialize_proteins_from_csv
from AFprep_func.uniprot_fetch import fetch_uniprot_info
from AFprep_func.classes import Domain
from AFprep_func.alphafold_db_domain_identification import read_afdb_json
from AFprep_func.alphafold_db_domain_identification import find_domains_from_pae
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

#adding domains from UniProt
for protein in proteins:
    data = fetch_uniprot_info(protein.accession_id)
    for feature in data['features']:
        if feature['type'] not in ['CHAIN', 'MUTAGEN', 'CONFLICT']:
            protein.add_domain(Domain(feature['description'], feature['begin'], feature['end'], feature['type']))

#adding user specified domains
df = pd.read_csv(CSV_PATH)
# Check if the "domains" column exists
if 'domains' in df.columns:
    # Convert the string representation of list of tuples in the "domains" column to actual lists of tuples
    df['domains'] = df['domains'].apply(lambda x: literal_eval(x) if pd.notnull(x) else [])

    # Create a name to protein mapping for easy lookup
    name_to_protein = {protein.name: protein for protein in proteins}

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        protein_name = row['name']
        # Check if this protein name exists in our proteins list
        if protein_name in name_to_protein:
            protein = name_to_protein[protein_name]
            # Retrieve domain data for this protein
            domain_data = row['domains']
            next_domain_num = 1  # Resetting or initializing domain numbering for each protein

            # Add domains to the protein
            for domain_start, domain_end in domain_data:
                protein.add_domain(Domain(f"manual_D{next_domain_num}", domain_start, domain_end, "manually_defined"))
                next_domain_num += 1
            print(next_domain_num, "manually specified domains added for ", protein.name)
else:
    print("The 'domains' column does not exist in the CSV file.")

#checking domains added
for protein in proteins:
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

#adding fragments to original csv file (both indices and sequences)
df = pd.read_csv(CSV_PATH)

fragments_dict = {protein.name: protein.fragments for protein in proteins}
df['fragments'] = df['name'].map(fragments_dict)

def get_fragment_list_sequence(sequence, fragment_list):
    fragments_sequence = [sequence[start:end] for start, end in fragment_list]
    return fragments_sequence
# Create a dictionary mapping protein names to their sequence fragments
fragments_sequence_dict = {protein.name: get_fragment_list_sequence(protein.sequence, protein.fragment_list) for protein in proteins}
df['fragments_sequence'] = df['name'].map(fragments_sequence_dict)

df.to_csv(CSV_PATH, index=False)
