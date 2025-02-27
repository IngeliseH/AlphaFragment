"""
This module provides a function to compile a list of Domain objects for a given
protein, using domain data from UniProt, AlphaFold structure predictions, and/or
manually specified domains.

Functions:
  - compile_domains: Compiles a list of Domain objects for a given protein

Dependencies:
  - .uniprot_fetch.find_uniprot_domains: For fetching domain data from
    UniProt.
  - .alphafold_db_domain_identification.read_afdb_json: For reading
    AlphaFold predicted structure data from the AFDB.
  - .alphafold_db_domain_identification.find_domains_from_pae: For
    identifying domains from AlphaFold predicted structures.
  - .process_proteins_csv.find_user_specified_domains: For adding
    manually specified domains from the input dataframe.
"""

def compile_domains(protein, uniprot=True, alphafold=True, manual=True,
                    protein_data=None, pae_method="definite", pae_custom_params=None):
    """
    Compiles a list of Domain objects for a given protein, using domain data
    from UniProt, AlphaFold structure predictions, and/or manually specified
    domains.

    Parameters:
      - proteins (list of Proteins): A list of Protein objects to be fragmented.
      - uniprot (bool, optional): Whether to fetch domain data from UniProt.
        Defaults to True.
      - alphafold (bool, optional): Whether to use domains identified from
        AlphaFold sructure predictions for the specified accession code, if
        available. Defaults to True.
      - manual (bool, optional): Whether to incorporate manually specified
        domains from the input CSV. Defaults to True.
      - protein_data (pandas.DataFrame): The dataframe containing protein data,
        including protein names and any manually specified domains. Required if
        manual is True. Defaults to None.
      - pae_method (str, optional): The method used to identify domains from
        AlphaFold predicted structures. Options are:
          - 'cautious' - Groups residues into domains with moderate PAE
            thresholds, aiming to balance sensitivity and specificity.
          - 'definite' - Only groups residues into domains if there is very high
            confidence in their relative positions. Likely to produce smaller
            domains.
          - 'custom' - Allows specification of custom PAE thresholds via the
            custom_params argument.
        Defaults to "definite".
      - pae_custom_params (dict, optional): Custom PAE thresholds for the AlphaFold
        domain identification method. Only use with pae_method 'custom". Should
        be a dictionary with the following keys:
            - 'res_dist_cutoff' (int) - The residue distance threshold to
            differentiate between close and further residue evaluations.
            Meaningless if close_pae_val == further_pae_val, or if set below 4
            as residues closer than this are ignored anyway due to high base
            confidence in their relative positions. Set to 10 for cautious
            grouping method, and irrelavant for definite grouping method as close
            and further pae_vals are equal.
          - 'close_pae_val' (int) - The PAE threshold which residue pairs must
            fall below to be considered within the same domain, if the distance
            between them is between 4 and res_dist_cutoff. Set to 4 for cautious
            grouping method and 2 for definite grouping method.
          - 'further_pae_val' (int) - The PAE threshold which residue pairs must
            fall below to be considered within the same domain, if the distance
            between them is greater than res_dist_cutoff. Set to 11 for close
            grouping method and 2 for definite grouping method.

    Returns:
      - List of Domain objects: A list of Domain objects, each representing a
        domain identified in the protein, with 'type' indicating the origin of
        the domain data ('UniProt', 'AlphaFold', 'manually_defined').
    
    Raises:
      - TypeError: If the protein argument is not an instance of the Protein class.

    Note:
      - Domain indexing is 0-based, so the start and end positions of the domains
        will be 1 less than standard residue positions. Manually provided domains
        are expected to be given in 1-based indexing.
    """
    from .classes import Protein
    from .uniprot_fetch import find_uniprot_domains
    from .alphafold_db_domain_identification import read_afdb_json, find_domains_from_pae
    from .process_proteins_csv import find_user_specified_domains

    #check that protein is instance of Protein
    if not isinstance(protein, Protein):
        raise TypeError("The 'protein' argument must be an instance of the Protein class.")
    if manual and protein_data is None:
        print("Cannot add manually specified domains - the 'df' argument is required.")
        manual = False
    elif manual and 'domains' not in protein_data.columns:
        print("Cannot add manually specified domains - the 'domains' column does not exist in the dataframe.")
        manual = False
    elif manual and 'name' not in protein_data.columns:
        print("Cannot add manually specified domains - the 'name' column does not exist in the dataframe, so domains cannot be matched with proteins.")
        manual = False

    domains = []

    # Identify and add domains from AlphaFold structure predictions
    if alphafold:
        protein_pae = read_afdb_json(protein.accession_id)
        if protein_pae:
            alphafold_domains = find_domains_from_pae(protein_pae, method=pae_method,
                                                      custom_params=pae_custom_params) or []
            if alphafold_domains:
                print(f"{len(alphafold_domains)} domains found in AlphaFold structure for "
                      f"{protein.name}: {alphafold_domains}")
                domains.extend(alphafold_domains)
            else:
                print(f"No domains found in AlphaFold structure for protein {protein.name}.")

    # Find and add domains from UniProt
    if uniprot:
        uniprot_domains = find_uniprot_domains(protein) or []
        domains.extend(uniprot_domains)

    # Add manually specified domains
    if manual:
        manual_domains = find_user_specified_domains(protein.name, protein_data) or []
        #print(f"{len(manual_domains)} domains found in csv for {protein.name}: {manual_domains}")
        domains.extend(manual_domains)

    return domains
