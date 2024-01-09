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
    target_id_data = target_id_data.append(data[col].str.split('; ').explode())

# Combine the "sample" column with the target ID columns
target_id_data['sample'] = samples.explode()

# Rename the column containing values from "_id" columns to "target_id"
target_id_data = target_id_data.rename(columns={0: 'target_id'})

# Create a new DataFrame for counting occurrences using crosstab
table = pd.crosstab(index=target_id_data['target_id'], columns=target_id_data['sample']).fillna(0)

# Save the resulting table to a TSV file
table.to_csv("target_id_counts.tsv", sep='\t')
