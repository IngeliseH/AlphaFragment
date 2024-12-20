===========
Usage Guide
===========

This section provides a quick start guide to using the AlphaFragment package.
This workflow demonstrates how to process protein data from a CSV file, identify
domains, fragment proteins, and visualize the results. For a more thorough guide
with explanations of how the package works and different options available, refer
to `AlphaFragment usage guide.ipynb <AlphaFragment_usage_guide.ipynb>`_.

Workflow
----------------------------

1. **Setting Up Your Python Script**

   Begin by setting up your Python script with the necessary paths for input and output:

   .. code-block:: python

       input_csv_path = "input.csv"
       output_csv_path = "output.csv"
       image_save_location = "folder_path"

2. **Initializing Protein Data from CSV**

   Import protein data from a CSV file. This step converts the CSV into a list
   of Protein objects and a DataFrame containing all data in the input file.

   .. code-block:: python

       from alphafragment import initialize_proteins_from_csv

       proteins, df = initialize_proteins_from_csv(input_csv_path)

3. **Protein Domain Identification and Fragmentation**

   For each protein in the dataset, identify domains and fragment the protein accordingly:

   .. code-block:: python

       from alphafragment import compile_domains, fragment_protein

       for protein in proteins:
           domains = compile_domains(protein, protein_data=df)
           for domain in domains:
               protein.add_domain(domain)

           fragments = fragment_protein(protein)
           for fragment in fragments:
               (protein.add_fragment(start, end) for start, end in fragment)

4. **Visualization of Fragmentation**

   Create a graphic that illustrates the domain locations and fragmentation results:

   .. code-block:: python
    
       from alphafragment import plot_fragmentation_output

       for protein in proteins:
           plot_fragmentation_output(protein, fragments, image_save_location)

5. **Updating and Saving Output Data**

   Update the DataFrame with protein information and save the updated data to a CSV file:

   .. code-block:: python

       from alphafragment import update_csv_with_fragments

       update_csv_with_fragments(df, output_csv_path, proteins)

6. **Output Generation**

   Generate FASTA files and AlphaPulldown input files for further analysis:

   .. code-block:: python

       from alphafragment import output_fastas, output_pulldown

       output_fastas(proteins)
       output_pulldown(proteins)
