import pandas as pd

# Read the data from the combined_annotations.tsv file
data = pd.read_csv("combined_annotations.tsv", sep='\t')

# Identify columns containing database IDs by searching for names ending with "_id"
id_columns = [col for col in data.columns if col.endswith("_id") and col != "query_id"]

# Split the "query_id" and "sample" columns by the delimiter "; " and create separate rows for each element
data['query_id'] = data['query_id'].str.split('; ')
data['sample'] = data['sample'].str.split('; ')

# Explode the lists to create separate rows for each element
data = data.explode('query_id')
data = data.explode('sample')

# Group by "query_id" and "sample," then count the occurrences for each database ID column
table = data.groupby(['query_id', 'sample'])[id_columns].count().fillna(0)

# Save the resulting table to a TSV file
table.to_csv("target_id_counts.tsv", sep='\t')
