import pandas as pd
import sys
import os

# Input arguments
annotations_file = sys.argv[1]
tree_option = sys.argv[2]  # Either "nar_nxr" or "amoA_pmoA"
output_file = sys.argv[3]

# Assuming the script is run from the Nextflow work directory
search_terms_file = f"{tree_option}/{tree_option}_search_terms.txt"

# Load search terms
with open(search_terms_file, 'r') as file:
    search_terms = [line.strip() for line in file]

# Load annotations
annotations = pd.read_csv(annotations_file, sep='\t')

# Filter annotations based on search terms and EC numbers
def filter_row(row):
    for term in search_terms:
        if term in row.to_string():
            return True
    return False

filtered_annotations = annotations[annotations.apply(filter_row, axis=1)]

# Extract query_ids
query_ids = filtered_annotations['query_id'].unique()

# Save query_ids to file
with open(output_file, 'w') as f:
    for query_id in query_ids:
        f.write(f"{query_id}\n")

print(f"Extracted {len(query_ids)} unique query IDs.")
