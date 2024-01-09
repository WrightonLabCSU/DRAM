import pandas as pd

# Read the data from the combined_annotations.tsv file
data = pd.read_csv("combined_annotations.tsv", sep='\t')

# Identify columns containing database IDs by searching for names ending with "_id"
id_columns = [col for col in data.columns if col.endswith("_id") and col != "query_id"]

# Extract "sample" column
samples = data['sample'].str.split('; ')

# Create a new DataFrame for target IDs and samples
target_id_data = pd.DataFrame()

# For each column ending with "_id", explode the list and create separate rows
for col in id_columns:
    target_id_data[col] = data[col].str.split('; ').explode()

# Combine the "sample" column with the target ID columns
target_id_data['sample'] = samples.explode()

# Create a new DataFrame for counting occurrences
occurrences_data = target_id_data.groupby(['sample'] + id_columns).size().reset_index(name='occurrences')

# Use the first column from target_id_data as the target_id values
occurrences_data['target_id'] = target_id_data.iloc[:, 0]

# Pivot the table to have samples as columns and target_ids as rows with the count of occurrences
table = occurrences_data.pivot(index='target_id', columns='sample', values='occurrences').fillna(0)

# Save the resulting table to a TSV file
table.to_csv("target_id_counts.tsv", sep='\t')
