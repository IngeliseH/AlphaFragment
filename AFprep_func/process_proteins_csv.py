"""
This module provides functionality for initializing protein objects from CSV
files, identifying domains in a CSV file associated with a particular protein
object, and updating the CSV with protein fragments information.

Functions:
  - initialize_proteins_from_csv: Initializes Protein objects from data in a CSV
    file.
  - add_user_specified_domains: Identifies domains listed in a CSV file for a
    corresponding Protein object.
  - update_csv_with_fragments: Updates a CSV file with information about protein
    fragments.

Dependencies:
  - pandas: For reading and processing the CSV file.
  - ast.literal_eval: For safely evaluating string literals containing Python
    expressions.
  - AFprep_func.classes.Protein: The Protein class for representing protein data.
  - AFprep_func.uniprot_fetch.fetch_uniprot_info: For fetching protein data from
    UniProt.
"""

import pandas as pd
from ast import literal_eval
# Import the Protein class
from AFprep_func.classes import Protein
# Import function to fetch sequence data
from AFprep_func.uniprot_fetch import fetch_uniprot_info

def initialize_proteins_from_csv(csv_path):
    """
    Reads a CSV file with headers for protein names and accession IDs, and
    initializes a list of Protein objects. Fetches protein sequences from
    UniProt - if fetching fails, can use manually provided sequence (eg if
    protein is known to not be in uniprot, provide a sequence in a column
    'sequence' and no accession_id - note that this will mean no domains can be
    found from alphafold or uniprot and these must be manually specified if
    needed). Does not initialize a Protein object for if both fetching and
    manual provision fail.

    Parameters:
      - csv_path (str): The path to the CSV file.

    Returns:
      - tuple of (list of Protein, list of str): 
          - The first element is a list of initialized Protein objects for which
            sequences were successfully fetched.
          - The second element is a list of protein names for which sequences
            could not be retrieved or other errors occurred during fetching.

    Errors and Exceptions:
      - Raises ValueError if required columns ('name', 'accession_id') are not
        found.
      - Reports when sequences cannot be fetched (e.g., invalid accession_id,
      network issues) and when a manually provided sequence is
        used instead.
      - Prints an error message if neither fetching nor manual sequence
        provision is successful and the protein cannot be initialised.
      
    Notes:
      - Sequences entered in a 'sequence' column will only be used if no
        sequence data is found in UniProt.
    """
    # Read the CSV file
    df = pd.read_csv(csv_path)

    # Normalize column names
    df.columns = [col.strip().lower().replace(" ", "_") for col in df.columns]

    # Check for required columns
    required_columns = ['name', 'accession_id']
    if not all(col in df.columns for col in required_columns):
        missing_cols = [col for col in required_columns if col not in df.columns]
        raise ValueError(f"Missing required columns: {', '.join(missing_cols)}")

    proteins = []  # List to store successfully initialized Protein objects
    proteins_with_errors = []  # List to store protein names with no UniProt data available

    # Iterate through each row in the DataFrame
    for _, row in df.iterrows():
        protein_name = row.get('name', '').strip()
        # Ensure each entry has a name; skip entries without a name
        if not protein_name:
            print("Missing protein name; skipping entry.")
            continue
        accession_id = row.get('accession_id', '')
        accession_id = '' if pd.isna(accession_id) else accession_id.strip()
        manual_sequence = row.get('sequence', '')
        manual_sequence = '' if pd.isna(manual_sequence) else manual_sequence.strip()

        sequence = ''
        # Attempt to fetch sequence
        fetch_result = fetch_uniprot_info(accession_id)
        if fetch_result and 'sequence' in fetch_result:
            sequence = fetch_result['sequence']
        elif manual_sequence:  # Use manually provided sequence if fetch fails
            print(f"UniProt fetch failed for {protein_name}; using manually provided sequence.")
            sequence = manual_sequence
        else:
            print(f"Unable to obtain a sequence for {protein_name} from UniProt or manual provision.")
            proteins_with_errors.append(protein_name)
            continue  # Skip this entry if no sequence available

        # Initialize Protein object if sequence is available
        proteins.append(Protein(name=protein_name, accession_id=accession_id, sequence=sequence))

    return proteins, proteins_with_errors

def add_user_specified_domains(protein, df):
    """
    Finds domains listed in a csv file for a particular protein.

    Parameters:
      - protein (Protein): Protein object to add domains to.
      - df (dataframe): DataFrame with protein names and user-defined domains.
    
    Notes:
      - Domains are identified by start and end positions within the protein sequence.
    """
    protein_name = protein.name
    manual_domains = []
    if protein_name in df['name'].values:
        domain_data = df.loc[df['name'] == protein_name, 'domains'].iloc[0]
        if pd.notnull(domain_data):
            domain_data = literal_eval(domain_data)
            for domain_start, domain_end in domain_data:
                next_domain_num = len(protein.domain_list) + 1
                manual_domains.append(Domain(f"manual_D{next_domain_num}", domain_start, domain_end, "manually_defined"))
            print(f"{len(domain_data)} manually specified domains found for {protein.name}: {manual_domains}")
    else:
        print(f"No user-specified domains found for protein {protein.name}.")
    return manual_domains

def update_csv_with_fragments(df, output_csv, proteins):
    """
    Reads a dataframe and outputs a csv file containing the original data, plus
    data for protein fragment indices and sequences, and domains identifed in
    the protein.

    Parameters:
      - df (dataframe): DataFrame containing initial information, likely from an
        input csv.
      - output_csv (file path ending .csv): Path to save the final csv file.
      - proteins (list of Protein objects): Protein objects containing
        information on Domain and fragment locations to be added to the output
        data.
    
    Notes:
      - Reorders columns for consistency
    """
    # Prepare dictionaries for mapping protein names to their updated attributes
    sequences_dict = {protein.name: protein.sequence for protein in proteins}
    domains_dict = {
        protein.name: ', '.join([f"{domain.num}: {domain.start}-{domain.end}" for domain in protein.domain_list])
        if protein.domain_list else ''
        for protein in proteins
    }
    print(domains_dict)
    fragments_dict = {protein.name: [(start, end) for start, end in protein.fragment_list] for protein in proteins}
    fragments_sequence_dict = {
        protein.name: [protein.sequence[start:end] for start, end in protein.fragment_list]
        for protein in proteins
    }

    # Update the DataFrame with the new sequences, domains, and fragment information
    df['sequence'] = df['name'].map(sequences_dict)
    df['domains'] = df['name'].map(domains_dict)
    df['fragment_indices'] = df['name'].map(fragments_dict)
    df['fragment_sequences'] = df['name'].map(fragments_sequence_dict)

    # Ensure the desired column order, adding other columns dynamically
    desired_columns = ['name', 'accession_id', 'sequence', 'domains', 'fragment_indices', 'fragment_sequences']
    existing_columns = df.columns.tolist()
    final_columns_order = [col for col in desired_columns if col in existing_columns] + [col for col in existing_columns if col not in desired_columns]
    # Reorder the DataFrame columns
    df = df.reindex(columns=final_columns_order)

    # Save the updated DataFrame to the new CSV file
    df.to_csv(output_csv, index=False)

    print(f"CSV file has been updated and saved to {output_csv}")
