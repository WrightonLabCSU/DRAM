import pandas as pd
import os
import argparse

def distill_summary(combined_annotations_path, distill_sheets_file, output_path):
    # Read the combined_annotations file
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')

    # Read the distill_sheets file and split paths
    with open(distill_sheets_file, 'r') as file:
        distill_sheets = [path.strip(',') for path in file.read().strip().split(',')]

    # Initialize an empty DataFrame to store the distill summary
    distill_summary_df = pd.DataFrame()

    # Process each distill sheet
    for distill_sheet in distill_sheets:
        # Print the path to the distill sheet for debugging
        print(f"Processing distill sheet: {distill_sheet}")

        # Read the distill sheet
        distill_df = pd.read_csv(distill_sheet, sep='\t')

        # Print the column names of the distill sheet for debugging
        print(f"Column names of distill sheet: {distill_df.columns}")

        # Identify the common gene_id column between combined_annotations_df and distill_df
        common_gene_id_columns = set(combined_annotations_df.columns) & set(distill_df.columns)

        # Exclude the "query_id" column from the common_gene_id_columns
        common_gene_id_columns = [col for col in common_gene_id_columns if col != "query_id"]

        # If gene_id is not found in the columns, try to identify by checking for a column ending with "_id"
        if not common_gene_id_columns:
            potential_gene_id_columns = [column for column in distill_df.columns if column.endswith('_id') and column != "query_id"]
            if potential_gene_id_columns:
                common_gene_id_columns = potential_gene_id_columns

        # If still not found, raise an error
        if not common_gene_id_columns:
            raise ValueError("No common gene_id column found between distill sheet and combined annotations.")

        # Merge the distill sheet with the combined_annotations using the common_gene_id_columns
        merged_df = pd.merge(combined_annotations_df, distill_df, left_on=common_gene_id_columns, right_on=common_gene_id_columns, how='inner')

        # Append the merged DataFrame to the distill summary DataFrame
        distill_summary_df = pd.concat([distill_summary_df, merged_df])

    # Save the distill summary to the specified output path
    distill_summary_df.to_csv(output_path, sep='\t', index=False, columns=['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory'])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--distill_sheets_file', required=True, help='Path to the distill_sheets file.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')

    args = parser.parse_args()

    distill_summary(args.combined_annotations, args.distill_sheets_file, args.output)
