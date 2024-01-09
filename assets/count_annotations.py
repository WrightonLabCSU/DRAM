import pandas as pd

# Read the data from the combined_annotations.tsv file
data = pd.read_csv("combined_annotations.tsv", sep='\t')

# Identify columns containing database IDs by searching for names ending with "_id"
id_columns = [col for col in data.columns if col.endswith("_id")]

# Extract "sample" column
samples = data['sample'].str.split('; ')

# Explode the "sample" list to create separate rows for each element
data = data.explode('sample')

# Create a new DataFrame for target IDs and samples
target_id_data = pd.DataFrame()

# For each column ending with "_id", explode the list and create separate rows
for col in id_columns:
    target_id_data[col] = data[col].str.split('; ').explode()

# Combine the "sample" column with the target ID columns
target_id_data['sample'] = samples.explode()

# Group by target IDs and samples, then count the occurrences
table = target_id_data.groupby(['sample'] + id_columns).size().unstack(fill_value=0)

# Save the resulting table to a TSV file
table.to_csv("target_id_counts.tsv", sep='\t')
