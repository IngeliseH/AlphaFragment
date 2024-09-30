from alphafragment.process_proteins_csv import initialize_proteins_from_csv, update_csv_with_fragments
from alphafragment.fragment_protein import fragment_protein
from alphafragment.domain_compilation import compile_domains
from alphafragment.plot_fragments import plot_fragmentation_output
from alphafragment.fragment_file_creation import output_fastas, output_pulldown

# specify input and output locations
input_csv_path = "input.csv"
output_csv_path = "output.csv"
image_save_location = "image_folder_path"
output_fastas_path = "fasta_folder_path"
output_pulldown_path = "pulldown_folder_path"
output_pulldown_fasta = "pulldown_fasta_path"

# Import protein data from a CSV file into a list of Protein objects and a DataFrame of all data in the CSV
proteins, df = initialize_proteins_from_csv(input_csv_path)

# Iterate through proteins to identify domains, fragment, and plot results
for protein in proteins:
    # Identify domains
    # Specify domain data sources to use as arguments - by default all are used
    domains = compile_domains(protein, protein_data=df)
    
    # Add domains to protein object
    for domain in domains:
        protein.add_domain(domain)

    # Fragment the protein with specified parameters
    fragments = fragment_protein(protein, length={'min': 150, 'ideal': 384, 'max': 400}, overlap={'min': 20, 'ideal': 30, 'max': 40})
    
    # Add fragments to protein object
    for fragment in fragments:
        protein.add_fragment(fragment)

    # Create graphic of domain locations and fragmentation results
    plot_fragmentation_output(protein, fragments, image_save_location, label=['UniProt', 'manually_defined'])

# Update DataFrame with protein information and save to CSV
update_csv_with_fragments(df, output_csv_path, proteins)

# Create FASTA files with fragment combinations
output_fastas(proteins, output_fastas_path)

# Alternatively, create AlphaPulldown input file with fragment combinations
output_pulldown(proteins, output_pulldown_path, fasta_name=output_pulldown_fasta)
