from AFprep_func.process_proteins_csv import initialize_proteins_from_csv, update_csv_with_fragments
from AFprep_func.fragment_protein import fragment_protein
from AFprep_func.compile_domains import compile_domains
from AFprep_func.plot_fragments import plot_fragmentation_output
from AFprep_func.fragment_file_creation import output_fastas, output_pulldown

input_csv_path = "input.csv"
output_csv_path = "output.csv"
image_save_location = "folder_path"

# import protein data from a csv file into a list of Protein objects and a dataframe of all data in the csv
proteins, df = initialize_proteins_from_csv(input_csv_path)

for protein in proteins:
    # identify domains
    domains = compile_domains(protein, protein_data=df)
    for domain in domains:
        protein.add_domain(domain)

    # fragment the protein
    fragments = fragment_protein(protein)
    for fragment in fragments:
        (protein.add_fragment(start, end) for start, end in fragment)

    # create graphic of domain locations and fragmentation results
    plot_fragmentation_output(protein, fragments, image_save_location, label=['UniProt', 'manually_defined'])

# update dataframe with protein information and save to csv
update_csv_with_fragments(df, output_csv_path, proteins)
# create fasta files with fragment combos
output_fastas(proteins)
# alternatively, create AlphaPulldown input file with fragment combos
output_pulldown(proteins)
