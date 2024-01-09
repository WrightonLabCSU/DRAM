import pandas as pd

# Read the data from the combined_annotations.tsv file
data = pd.read_csv("combined_annotations.tsv", sep='\t')

# Identify columns containing database IDs by searching for names ending with "_id"
id_columns = [col for col in data.columns if col.endswith("_id") and col != "query_id"]

# Extract "sample" column
samples = data['sample'].str.split('; ')

# Create an empty DataFrame to store the target_id data
target_id_data = pd.DataFrame()

# Process the input annotation files
for i in range(0, len(annotation_files), 2):
    # Read the data from the annotation file
    data = pd.read_csv(annotation_files[i + 1], sep='\t')

    # Identify columns containing database IDs by searching for names ending with "_id"
    id_columns = [col for col in data.columns if col.endswith("_id") and col != "query_id"]

    # Create a DataFrame with 'sample' and the current '_id' column exploded
    exploded_data = pd.concat([data[['sample', col]].explode(col) for col in id_columns])

    # Append the exploded data to the target_id_data DataFrame
    target_id_data = pd.concat([target_id_data, exploded_data])

# Create a new DataFrame for counting occurrences using groupby
table = target_id_data.groupby(['sample', 'target_id']).size().unstack(fill_value=0)

# Save the resulting table to a TSV file
table.reset_index().to_csv("target_id_counts.tsv", sep='\t', index=False)
