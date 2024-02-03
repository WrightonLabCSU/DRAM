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

    # Flag to check if "potential_amg" exists in any input sheet
    include_potential_amg = False

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

        # Check if the "potential_amg" column exists in the distill sheet
        if 'potential_amg' in distill_df.columns:
            include_potential_amg = True

        # Initialize an empty DataFrame to store the merged data for the current distill sheet
        merged_data_for_current_gene_id = pd.DataFrame()

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
            merged_data_for_current_gene_id = pd.concat([merged_data_for_current_gene_id, merged_df])

        # Merge with target_id_counts based on 'gene_id' and 'target_id'
        merged_data_for_current_gene_id = pd.merge(merged_data_for_current_gene_id, target_id_counts_df, left_on=['gene_id'], right_on=['target_id'], how='left')

        # Append the merged data for the current distill sheet to the overall distill summary DataFrame
        distill_summary_df = pd.concat([distill_summary_df, merged_data_for_current_gene_id])

    # Append "potential_amg" column to columns_to_output if it exists
    columns_to_output = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']
    if include_potential_amg:
        columns_to_output.append('potential_amg')

    # Extract the sample names from the target_id_counts columns (excluding non-numeric columns)
    sample_columns = target_id_counts_df.columns[target_id_counts_df.dtypes == 'int64']
    sample_names = sample_columns.tolist()

    # Append the sample-named columns to the columns_to_output
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
