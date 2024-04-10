"""
This module encapsulates the process of protein data fragmentation in a
domain-aware manner.

Functions:
  - fragment_proteins_from_csv: The main function that orchestrates the reading
    of protein data from a CSV file, the addition of domain data from various
    sources, fragmentation of the protein whtout breaking the identified domains,
    and the output of the resulting fragments to a new CSV file along with
    optional graphics showing fragment locations.

Dependencies:
  - pandas: For CSV file manipulation.
  - AFprep_func submodules: A collection of functions from the AFprep_func package for
    domain identification, protein fragmentation, and output visualization.
"""
import pandas as pd
from AFprep_func.uniprot_fetch import find_uniprot_domains
from AFprep_func.alphafold_db_domain_identification import read_afdb_json, find_domains_from_pae
from AFprep_func.fragment_protein import fragment_protein
from AFprep_func.plot_fragments import plot_fragmentation_output
from AFprep_func.process_proteins_csv import initialize_proteins_from_csv, add_user_specified_domains, update_csv_with_fragments

def fragment_proteins_from_csv(input_csv_path, output_csv_path, uniprot=True,
                               alphafold=True, manual=True, image_save_location=None):
    """
    Orchestrates the domain aware fragmentation of proteins listed in an input
    CSV file, using domain data from a specified combination of UniProt,
    AlphaFold, and manual specifications. Outputs the fragmented proteins along
    with any additional data in the original csv to a new CSV file. Optionally,
    saves visualisations of the fragmentation.

    Parameters:
      - input_csv_path (str): Path to the input CSV file containing protein data.
      - output_csv_path (str): Path where the output CSV file will be saved,
        including the original data plus domain and fragmentation information.
      - uniprot (bool, optional): Whether to fetch domain data from UniProt.
        Defaults to True.
      - alphafold (bool, optional): Whether to use domains identified from
        AlphaFold sructure predictions for the specified accession code, if
        available. Defaults to True.
      - manual (bool, optional): Whether to incorporate manually specified
        domains from the input CSV. Defaults to True.
      - image_save_location (str, optional): Directory path where graphics of
        protein fragmentation will be saved. If None, images are produced but
        not saved. Defaults to None.

    Returns:
      None. The function produces an output CSV file at `output_csv_path` and,
      optionally, graphics of the domain locations and fragmentation results.

    Notes:
      - For proteins with no UniProt data available, the function will use a
        sequence specified in the CSV file and any domains specified manually in
        the CSV. For protein for which UniProt data is available, and sequence
        specified in the CSV file will be overwritten by the UniProt sequence.
    """
    # read in data from csv
    proteins, errors = initialize_proteins_from_csv(input_csv_path)
    print(f"Successfully initialized proteins: {[protein.name for protein in proteins]}")
    print(f"Proteins with errors or no data available: {errors}")

    # add domains
    df = pd.read_csv(input_csv_path)

    if manual and 'domains' not in df.columns:
        print("Cannot add manually specified domains - the 'domains' column does not exist in the CSV file.")
        manual = False

    for protein in proteins:
        domains = []
        if uniprot:
            uniprot_domains = find_uniprot_domains(protein) or []
            domains.extend(uniprot_domains)
        if alphafold:
            protein_pae = read_afdb_json(protein.accession_id)
            if protein_pae:
                alphafold_domains = find_domains_from_pae(protein_pae) or []
                if alphafold_domains:
                    print(f"{len(alphafold_domains)} domains found in AlphaFold structure for {protein.name}: {alphafold_domains}")
                    domains.extend(alphafold_domains)
                else:
                    print(f"No domains found in AlphaFold structure for protein {protein.name}.")
        if manual:
            manual_domains = add_user_specified_domains(protein, df) or []
            domains.extend(manual_domains)

        for domain in domains:
            protein.add_domain(domain)

        #fragment protein
        fragments = fragment_protein(protein)
        print("Cutpoints found for protein ", protein.name, ":", fragments)
        for start, end in fragments:
            protein.add_fragment(start, end)
        plot_fragmentation_output(protein, protein.fragment_list, image_save_location)
        print(protein.name, " is ", len(protein.sequence), " residues long and has ",
            len(protein.fragment_list), "fragments")

    # update dataframe with protein information and save to csv
    update_csv_with_fragments(df, output_csv_path, proteins)
