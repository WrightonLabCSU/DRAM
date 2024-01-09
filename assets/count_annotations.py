import argparse
import pandas as pd

def count_annotations(input_file, output_file):
    # Read the data from the combined_annotations.tsv file
    data = pd.read_csv(args.input_file, sep='\t')

    # Identify columns containing database IDs by searching for names ending with "_id"
    id_columns = [col for col in data.columns if col.endswith("_id") and col != "query_id"]

    # Create an empty DataFrame to store the combined data
    target_id_data = pd.DataFrame(columns=['sample', 'target_id'])

    # Process each ID column
    for col in id_columns:
        # Explode the lists to create separate rows for each element
        exploded_data = data[['sample', col]].explode(col)
        # Rename the columns for consistency
        exploded_data.columns = ['sample', 'target_id']
        # Append the exploded data to the target_id_data DataFrame
        target_id_data = pd.concat([target_id_data, exploded_data], ignore_index=True)

    # Group by "target_id" and "sample," then count the occurrences
    table = target_id_data.groupby(['target_id', 'sample']).size().unstack(fill_value=0)

    # Save the resulting table to a TSV file
    table.to_csv(args.output_file, sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count occurrences of target_ids for each sample.")
    parser.add_argument("input_file", help="Input file (combined_annotations.tsv)")
    parser.add_argument("output_file", help="Output file (target_id_counts.tsv)")

    args = parser.parse_args()

    count_annotations(args.input_file, args.output_file)
