"""
Test file for the end-to-end workflow of the package
"""
import os
import pandas as pd
from alphafragment.process_proteins_csv import initialize_proteins_from_csv, update_csv_with_fragments
from alphafragment.fragment_protein import fragment_protein
from alphafragment.domain_compilation import compile_domains
from alphafragment.plot_fragments import plot_fragmentation_output
from alphafragment.fragment_file_creation import output_fastas, output_pulldown


def test_end_to_end_workflow():
    """
    Test the end-to-end workflow of the package, from reading an input CSV to generating
    output files and visualizations.
    """
    # Setup
    input_csv_path = '../sample_data/sample_input.csv'
    output_csv_path = '../sample_data/sample_output.csv'
    output_fastas_path = '../sample_data/sample_output_fastas'
    output_pulldown_path = '../sample_data/sample_pulldown.txt'
    image_save_location = '../sample_data/sample_output_images'
    os.makedirs(image_save_location, exist_ok=True)

    # Prepare a mock input CSV
    df_mock = pd.DataFrame({
        'name': ['ProteinA', 'ProteinB'],
        'accession_id': ['A', 'B'],
        'sequence': ['DDDDDXDDDDDDXXXXXXXDDDDDDDDDDDDDDDDXXXXDDDDDDXXXXX', 'IIIIIIIIIIIIIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIII'],
        'domains': [[(1, 5), (7, 12), (20, 35), (40, 45)], [(15, 45)]]
    })
    df_mock.to_csv(input_csv_path, index=False)

    # Run the workflow
    proteins, df = initialize_proteins_from_csv(input_csv_path)
    for protein in proteins:
        domains = compile_domains(protein, protein_data=df)
        for domain in domains:
            protein.add_domain(domain)
        length = {'min': 10, 'ideal': 15, 'max': 20}
        overlap = {'min': 3, 'ideal': 5, 'max': 7}
        fragments = fragment_protein(protein, length, overlap)
        print(f"fragments = {fragments}")
        for fragment in fragments:
            protein.add_fragment(fragment)
        plot_fragmentation_output(protein, fragments, image_save_location, label=['UniProt', 'manually_defined'])

    update_csv_with_fragments(df, output_csv_path, proteins)
    output_fastas(proteins, output_fastas_path)
    output_pulldown(proteins, output_pulldown_path)

    # Assertions
    assert os.path.exists(output_csv_path), "Output CSV file was not created"
    assert os.path.exists(output_pulldown_path), "Output pulldown file was not created"
    assert len(os.listdir(image_save_location)) == len(proteins), "Not all protein images were generated"

    # Cleanup
    os.remove(input_csv_path)
    os.remove(output_csv_path)
    for file in os.listdir(image_save_location):
        os.remove(os.path.join(image_save_location, file))
    os.rmdir(image_save_location)
    os.remove(output_pulldown_path)
    for file in os.listdir(output_fastas_path):
        os.remove(os.path.join(output_fastas_path, file))
    os.rmdir(output_fastas_path)
