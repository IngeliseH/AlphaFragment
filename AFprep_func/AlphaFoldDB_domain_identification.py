#identifying domains from AlphaFold DB

#reading in list of accession codes in AlphaFold DB - should use to give exception if AlphaFoldDB does not have prediction for input id
curl http://ftp.ebi.ac.uk/pub/databases/alphafold/accession_ids.txt -o accession_ids.txt
import requests
url = "http://ftp.ebi.ac.uk/pub/databases/alphafold/accession_ids.txt"
response = requests.get(url)
print(response)

accession_contents = open('accession_ids.txt', 'r').readlines()
accession_dict = {}
for i in accession_contents:
    accession_dict[i.split(',')[0]] = i.split(',')[3]
accession_id = 'P14618'
alphafold_id = None
if accession_id in accession_dict:
    alphafold_id = accession_dict[accession_id]