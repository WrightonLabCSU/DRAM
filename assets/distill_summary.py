import pandas as pd
import os
import argparse
import glob
import logging

logging.basicConfig(level=logging.DEBUG)

def is_null_content(file_path):
    with open(file_path, 'r') as file:
        content = file.read().strip()
    return content == "NULL"

def distill_summary(combined_annotations_path, target_id_counts_df, output_path):
    # Read the combined_annotations file
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')

    # Identify the potential gene ID columns in combined_annotations_df
    potential_gene_id_columns = [col for col in combined_annotations_df.columns if col.endswith('_id') and col != "query_id"]

    # Get a list of all _distill_sheet.tsv files in the current working directory
    distill_sheets = glob.glob('*_distill_sheet.tsv')

    # Initialize an empty DataFrame to store the distill summary
    distill_summary_df = pd.DataFrame()

    # Initialize an empty list for sample names
    sample_names = []

    # Initialize a dictionary to store additional columns dynamically
    additional_columns = {}

    # Process each distill sheet
    for distill_sheet in distill_sheets:
        # Check if the distill sheet content is "NULL"
        if is_null_content(distill_sheet):
            logging.info(f"Skipping distill sheet '{distill_sheet}' as it contains 'NULL' content.")
            continue

        # Print the path to the distill sheet for debugging
        logging.info(f"Processing distill sheet: {distill_sheet}")

        # Read the distill sheet
        distill_df = pd.read_csv(distill_sheet, sep='\t')

        # Include "potential_amg" column in columns_to_output if it exists
        columns_to_output = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']

        # Check if additional columns exist in the distill sheet
        additional_cols = [col for col in distill_df.columns if col not in columns_to_output]
        
        # Process each potential gene ID column
        for common_gene_id_column in potential_gene_id_columns:
            # Merge the distill sheet with the combined_annotations using the current gene ID column
            merged_df = pd.merge(
                combined_annotations_df,
                distill_df,
                left_on=[common_gene_id_column],
                right_on=['gene_id'],
                how='inner'
            )

            # Append the merged DataFrame to the distill summary DataFrame for the current distill sheet
            distill_summary_df = pd.concat([distill_summary_df, merged_df])

            # Update the additional columns dictionary with matching gene IDs
            for col in additional_cols:
                if col not in additional_columns:
                    additional_columns[col] = {}
                additional_columns[col].update(merged_df.set_index('gene_id')[col].to_dict())

        # Merge with target_id_counts based on 'gene_id' and 'target_id'
        distill_summary_df = pd.merge(distill_summary_df, target_id_counts_df, left_on=['gene_id'], right_on=['target_id'], how='left')


    # Append additional columns and their values to columns_to_output
    for col, values in additional_columns.items():
        distill_summary_df[col] = distill_summary_df['gene_id'].map(values)
        columns_to_output.append(col)

    # Extract the sample names from the target_id_counts columns (excluding non-numeric columns)
    sample_columns = target_id_counts_df.columns[target_id_counts_df.dtypes == 'int64']
    sample_names = sample_columns.tolist()

    # Append the sample-named columns to columns_to_output
    columns_to_output += sample_names

    # Save the deduplicated distill summary to the specified output path
    deduplicated_df = distill_summary_df.drop_duplicates(subset=['gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory'])
    deduplicated_df.to_csv(output_path, sep='\t', index=False, columns=columns_to_output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--target_id_counts', required=True, help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')

    args = parser.parse_args()

    # Read the target_id_counts file
    target_id_counts_df = pd.read_csv(args.target_id_counts, sep='\t')

    distill_summary(args.combined_annotations, target_id_counts_df, args.output)
