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