"""
This module provides functionality for initializing protein objects from CSV
files, identifying domains in a CSV file associated with a particular protein
object, and updating the CSV with protein fragments information.

Functions:
  - initialize_proteins_from_csv: Initializes Protein objects from data in a CSV
    file.
  - find_user_specified_domains: Identifies domains listed in a CSV file for a
    corresponding Protein object.
  - update_csv_with_fragments: Updates a CSV file with information about protein
    fragments.

Dependencies:
  - pandas: For reading and processing the CSV file.
  - ast.literal_eval: For safely evaluating string literals containing Python
    expressions.
  - collections.Counter: For checking column names are not duplicated.
  - .classes.Protein: The Protein class for representing protein data.
  - .classes.Domain: The Domain class for representing protein domains.
  - .uniprot_fetch.fetch_uniprot_info: For fetching protein data from UniProt.
"""
from ast import literal_eval
from collections import Counter
import pandas as pd
from .classes import Protein, Domain
from .uniprot_fetch import fetch_uniprot_info

def initialize_proteins_from_csv(csv_path):
    """
    Reads a CSV file with columns for protein names and accession IDs, and
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
      - Column names must be 'name' and 'accession_id' for protein name and
        accession ID respectively. If csv includes columns for manually specified
        protein sequences or domains, these should be named 'sequence' and
        'domains' respectively. Capitalisation will be ignored. Other columns can
        also be included.
      - Sequences entered in a 'sequence' column will only be used if no
        sequence data is found in UniProt. If a sequence is found in UniProt, the
        'sequence' column will be ignored (and overwritten if the
        update_csv_with_fragments function is used).
    """
    # Read the CSV file
    df = pd.read_csv(csv_path)

    # Normalize column names
    df.columns = [col.strip().lower().replace(" ", "_") for col in df.columns]

    # Check for duplicate column names after normalization
    column_counter = Counter(df.columns)
    duplicates = [col for col, count in column_counter.items() if count > 1]
    if duplicates:
        raise ValueError("Duplicate column names detected after normalization: "
                         f"{', '.join(duplicates)}")

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

        # Initialize Protein object if sequence is available
        if sequence:
            proteins.append(Protein(name=protein_name, accession_id=accession_id,
                                    sequence=sequence))
        else:
            print(f"Unable to obtain a sequence for {protein_name} from UniProt "
                  f"or manual provision. Cannot initialize Protein object.")
            proteins_with_errors.append(protein_name)

    if proteins:
        print(f"Successfully initialized proteins: {[protein.name for protein in proteins]}")
    if proteins_with_errors:
        print(f"Proteins with errors or no data available: {proteins_with_errors}")
    return proteins, df

def find_user_specified_domains(protein_name, df):
    """
    Finds domains listed in a dataframe for a particular protein.

    Parameters:
      - protein_name (str): The name of the protein to find domains for.
      - df (dataframe): DataFrame with protein names in a column 'name' and
        user-defined domains in a column "domains", as a list of (start, end)
        tuples for each protein, using 1-based indexing.

    Returns:
      - list of Domain objects: A list of Domain objects representing the
        manually defined domains for the protein, using 0-based indexing.

    Raises:
      - TypeError: If the input DataFrame or domain_data has incorrect type.
      - ValueError: If the DataFrame is missing required columns or if the domain
        data cannot be parsed.

    Notes:
      - Domains are expected to be provided using 1-based indexing but will be
        processed and output with 0-based indexing.
    """
    def create_domain(identifier, start, end):
        """ Helper to create a Domain object, converting from 1-based to 0-based indexing."""
        if not (isinstance(start, int) and isinstance(end, int)):
            raise TypeError("Each domain must consist of two integers, but "
                            f"received: {start}, {end}")
        return Domain(identifier, start-1, end-1, "manually_defined")

    # Validate DataFrame input
    if not isinstance(df, pd.DataFrame):
        raise TypeError("The 'df' argument must be a pandas dataframe.")
    if 'name' not in df.columns or 'domains' not in df.columns:
        missing_cols = [col for col in ['name', 'domains'] if col not in df.columns]
        raise ValueError(f"Missing columns in dataframe: {', '.join(missing_cols)}")

    # Filter DataFrame for the specified protein
    protein_domains = df[df['name'] == protein_name]
    # Return empty list if protein is not found
    if protein_domains.empty:
        print(f"No user-specified domains found for protein {protein_name}.")
        return []

    # Return empty list if protein is found but has no associated domains
    domain_data = protein_domains['domains'].iloc[0]
    if not domain_data:
        print(f"No user-specified domains found for protein {protein_name}.")
        return []

    # Validate and parse domain data if it's a string representation of a list
    if isinstance(domain_data, str):
        try:
            domain_data = literal_eval(domain_data)
        except (SyntaxError, ValueError) as e:
            raise ValueError(f"Error parsing domain data: {str(e)}") from e

    domains = []
    # Create Domain objects from the domain data (list of tuples or dictionary)
    if isinstance(domain_data, list):
        domains = ([create_domain(f"manual_D{index + 1}", *domain) for index, domain in enumerate(domain_data)])
    if isinstance(domain_data, dict):
        domains = ([create_domain(domain_id, *positions) for domain_id, positions in domain_data.items()])
    
    if domains:
        print(f"{len(domains)} domains found in csv for {protein_name}: {domains}")
        return domains

    # Raise an error if the domain data is not a list or dictionary
    raise TypeError("Domain data must be a list of (start, end) tuples or a "
                    "dictionary of (start, end) tuples with identifiers.")


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

    Returns:
      - dataframe: A new DataFrame with updated data and reordered columns.
    
    Notes:
      - Reorders columns for consistency
      - Domains in output are referenced using 1-based indexing.
      - Fragments in output are referenced using 1-based indexing and inclusive
        of the start and end residues.
    """
    # Create a copy of the DataFrame to avoid modifying the original
    df_copy = df.copy()

    # Prepare dictionaries for mapping protein names to their updated attributes
    sequences_dict = {protein.name: protein.sequence for protein in proteins}
    domains_dict = {
        protein.name: ', '.join([f"{domain.id}: {domain.start+1}-{domain.end+1}"
                                 for domain in protein.domain_list])
        if protein.domain_list else ''
        for protein in proteins
    }
    fragments_dict = {
        protein.name: [(start + 1, end) for start, end in protein.fragment_list]
        for protein in proteins
    }
    fragments_sequence_dict = {
        protein.name: [protein.sequence[start:end] for start, end in protein.fragment_list]
        for protein in proteins
    }

    # Update the DataFrame copy with the new sequences, domains, and fragment information
    df_copy['sequence'] = df['name'].map(sequences_dict)
    df_copy['domains'] = df['name'].map(domains_dict)
    df_copy['fragment_indices'] = df['name'].map(fragments_dict)
    df_copy['fragment_sequences'] = df['name'].map(fragments_sequence_dict)

    # Define desired column order, adding other columns dynamically
    desired_columns = ['name', 'accession_id', 'sequence', 'domains',
                       'fragment_indices', 'fragment_sequences']
    existing_columns = df_copy.columns.tolist()
    additional_columns = [col for col in existing_columns if col not in desired_columns]
    final_columns_order = ([col for col in desired_columns if col in existing_columns]
                           + additional_columns)

    # Reorder the DataFrame columns
    new_df = df_copy.reindex(columns=final_columns_order)
    # Save the updated DataFrame to the new CSV file
    new_df.to_csv(output_csv, index=False)
    print(f"CSV file has been updated and saved to {output_csv}")

    # Return the newly organized DataFrame
    return new_df
