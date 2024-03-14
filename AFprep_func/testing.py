from protein_init import initialize_proteins_from_csv  # Importing the function to fetch sequence data

# Example usage
csv_path = 'AFprep_func/proteins1.csv'
proteins, errors = initialize_proteins_from_csv(csv_path)
print(f"Successfully initialized proteins: {[protein.name for protein in proteins]}")
print(f"Proteins with errors or no data available: {errors}")

###############
#Domain identification from AF
from AlphaFoldDB_domain_identification import read_AFDB_json
from AlphaFoldDB_domain_identification import find_domains_from_PAE


#testing using Ana2
Ana2_PAE = read_AFDB_json("Q9XZ31")
print(Ana2_PAE)

# testing using Ana2
Ana2_domains = find_domains_from_PAE(Ana2_PAE)
print(Ana2_domains)
print(len(Ana2_PAE))

# testing function by visualising output for Ana2
import pandas as pd
import matplotlib
matplotlib.use('Qt5Agg', force=True)  # Force the use of the Qt5Agg backend

# Now import pyplot after setting the backend
import matplotlib.pyplot as plt


# Convert the dictionary to a list of tuples [(key, item), ...] for easier DataFrame creation
data = [(key, item) for key, value_list in Ana2_domains.items() for item in value_list]
print(data)

# Create the DataFrame
df = pd.DataFrame(data, columns=['Domain', 'Item'])
print(df)

# Plotting
plt.figure(figsize=(10, 6))
# Assuming 'Item' should be treated as a continuous variable on the x-axis and 'Domain' as a categorical variable on the y-axis
# Convert domain to category type for better handling by matplotlib
df['Domain'] = pd.Categorical(df['Domain'])

# Scatter plot with 'Item' on x-axis and 'Domain' on y-axis
plt.scatter(df['Item'], df['Domain'])

plt.xlabel('Item')
plt.ylabel('Domain')
plt.title('Domain Distribution of Items')
plt.xlim(0, len(Ana2_PAE))  # Set x-axis limits based on the length of Ana2_PAE
plt.grid(True)

plt.show()


#######
#testing recursive fragmentation on mock protein
from recursive_fragment_finding import recursive_fragmentation
from classes import Domain
from AlphaFoldDB_domain_identification import read_AFDB_json
from AlphaFoldDB_domain_identification import find_domains_from_PAE


# Example call (values for first_res, last_res, min_len, max_len need to be defined)
first_res = 1
last_res = 1000  # Example end of sequence
domains = [
    Domain('D1', 1, 190, 'manual'),
    Domain('D2', 400, 430, 'manual'),
    Domain('D3', 435, 500, 'manual'),
    Domain('D4', 300, 350, 'manual'),
    Domain('D5', 565, 700, 'manual')
]
min_len = 100
max_len = 200
overlap = 10
max_overlap = 20
min_overlap = 0
result = recursive_fragmentation(first_res, first_res, last_res, min_len, max_len, overlap, max_overlap, min_overlap)

if result is not None:
    print("Cutpoints found:", result)
else:
    print("No valid fragmentation pattern found.")

###########
#testing recursive fragmentation on real proteins
from protein_init import initialize_proteins_from_csv  # Importing the function to fetch sequence data
# Example usage
csv_path = 'proteins2.csv'
proteins, errors = initialize_proteins_from_csv(csv_path)
print(f"Successfully initialized proteins: {[protein.name for protein in proteins]}")
print(f"Proteins with errors or no data available: {errors}")

for protein in proteins:
    protein_PAE = read_AFDB_json(protein.accession_id)
    domains = find_domains_from_PAE(protein_PAE)
    for domain in domains:
        protein.add_domain(domain)
    print(protein.name, " has ", len(protein.domain_list), "domains")

min_len = 100
max_len = 200
overlap = 10
max_overlap = 20
min_overlap = 0
for protein in proteins:
    first_res = 0
    last_res = len(protein.sequence)
    while fragments == None:
        fragments = recursive_fragmentation(first_res, first_res, last_res, min_len, max_len, overlap, max_overlap, min_overlap)
        if fragments == None:
            min_len+=10
            max_len+=10
    if fragments is not None:
        print("Cutpoints found for protein ", protein, ":", fragments)
    else:
        print("No valid fragmentation pattern found for protein ", protein, ":")
