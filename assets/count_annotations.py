import pandas as pd

# Read the data from the combined_annotations.tsv file
data = pd.read_csv("combined_annotations.tsv", sep='\t')

# Split the "target_id" and "sample" columns by the delimiter "; " and create separate rows for each element
data['target_id'] = data['target_id'].str.split('; ')
data['sample'] = data['sample'].str.split('; ')

# Explode the lists to create separate rows for each element
data = data.explode('target_id')
data = data.explode('sample')

# Group by "target_id" and "sample," then count the occurrences
table = data.groupby(['target_id', 'sample']).size().unstack(fill_value=0)

# Save the resulting table to a TSV file
table.to_csv("target_id_counts.tsv", sep='\t')
