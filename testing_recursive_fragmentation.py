#testing recursive fragmentation on real protein
from AFprep_func.classes import Protein
from AFprep_func.alphafold_db_domain_identification import read_afdb_json
from AFprep_func.alphafold_db_domain_identification import find_domains_from_pae
from AFprep_func.protein_init import initialize_proteins_from_csv  # Importing the function to fetch sequence data
from AFprep_func.fragment_protein import fragment_protein

csv_path = 'AFprep_func/proteins2.csv'
#csv_path = 'proteins2.csv'
proteins, errors = initialize_proteins_from_csv(csv_path)
print(f"Successfully initialized proteins: {[protein.name for protein in proteins]}")
print(f"Proteins with errors or no data available: {errors}")

for protein in proteins:
    print(protein)
proteins.append(Protein(name='Plp',
                        accession_id='Q400N1',
                        sequence='MQTAEVSKMFQGPKIPLRELQHERDTLQGRMEEQTQRISTLQNRLEEQRQRAEQLHRTGTSDLNTRVHELQGEVQNLYEQLAARDKQMANMRQQLQRSKEEITRLETEVEVRTQPDRSLVNKLQAEVQQKGAEIVKLKDKIRTEMINRLAIPDLMETMLADKNDEIDHLRDQLEAKEKELQASQQEASQISSPSGAAGKQEGSGGKLSARTLSDIGSITEFPEPDVERRAAMRSLTALQMSEGAGGFLHQTMETSKEAVANLTHKRTDDLSGFIVPYPVNTFEHPHYFQALGVTAQSTDGLTPGLVPRQINFSNLTEDSKLKTTSLLMHTPELPRPTTPPEIHQLRVKLSDLQTEKQRQQSELEKKLQDLQKELEQEKEKLSRQAQTLQSYEESEAKYRLRIENLESKVLETAAQAASDRENLRKELNCVSAAHEQCENAAAARKRELEKLNSEVKVKADQLHAALRRCADLELQVLTLERDLERLKNSDNSSKQYSVDEIAQQVEKELNYSAQLDSNILKAIESEEENNLDKKLQKGVQTEEETLPGTGNGTDDENFTGERELLNQLEALRAQLAVEREQCEAMSKELLGEKQHSQDIQEQDVIIIEAMRKRLETALDAEDELHKQLDQERERCERLQTQLTSLQRAESRRNSSLLLKSPGDSPRKSPRADFESELGDRLRSEIKLLVAQNERERERSADAQRSSERERQRYEKELQERVAYCERLKQEMEKLSRDKESAETELEHFNERLTLQASEIESLEARLVTLQEAETRRANTRTRQHQENVKLQAEIHELKSKLLAAEAARDCLDQKVTQLRFDVSRSGQREAKLAEALAQANDRLAHSTDDNVPAQFMQKMKEINALLAENTQENRQMAETVQFLVGERIALQKKCEELGGAGNTNVTELEERCRQLIGRYLRVESHRKALVYQKRYLKLTLEGYQASEQLALQNLRGGAEQPPQRNIKKKFKTVALAIIAIQRIKYIGRIWHTGKRIVSKSVFTITQQRSPGLNLNVAPPQSPLPGPNSNLPTNNNNLMGRLSYAPLSPPMVNFSTVQPIMLPSDFTLQAPTTSLQTTIANNNSENTLPSLARLDWPTMQKPKRAHARHH'))

for protein in proteins:
    protein_PAE = read_afdb_json(protein.accession_id)
    domains = find_domains_from_pae(protein_PAE)
    for domain in domains:
        print(domain)
        protein.add_domain(domain)
    print(protein.name, " has ", len(protein.domain_list), "domains")

    fragment_protein(protein)


