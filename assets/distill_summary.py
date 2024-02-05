import pandas as pd
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

    # Initialize a set to store unique additional columns dynamically
    additional_columns = set()

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

        # Process each potential gene ID column
        for common_gene_id_column in potential_gene_id_columns:
            # Merge the distill sheet with the combined_annotations using the current gene ID column
            merged_df = pd.merge(
                combined_annotations_df,
                distill_df,
                left_on=[common_gene_id_column],
                right_on=['gene_id'],
                how='inner',
                suffixes=('', '_dup')
            )

            # Remove the duplicate columns resulted from merging
            for col in merged_df.columns:
                if col.endswith('_dup'):
                    merged_df.drop(col, axis=1, inplace=True)

            # Append the merged DataFrame to the distill summary DataFrame for the current distill sheet
            distill_summary_df = pd.concat([distill_summary_df, merged_df], ignore_index=True)

        # Update the set of additional columns with those not in combined_annotations_df
        additional_columns.update(set(distill_df.columns) - set(combined_annotations_df.columns) - {'gene_id'})

    # Merge with target_id_counts based on 'gene_id' and 'target_id', including the additional columns
    distill_summary_df = pd.merge(distill_summary_df, target_id_counts_df, left_on=['gene_id'], right_on=['target_id'], how='left')

    # Define columns_to_output with the basic columns present
    basic_columns = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']

    # Deduplicate the basic columns
    distill_summary_df.drop_duplicates(subset=basic_columns, inplace=True)

    # Add the dynamically determined additional columns
    columns_to_output = basic_columns + list(additional_columns)

    # Add the dynamically determined bin columns from the target_id_counts_df
    bin_columns = [col for col in target_id_counts_df.columns if col.startswith('bin-')]
    columns_to_output.extend(bin_columns)

    # Save the deduplicated distill summary to the specified output path
    distill_summary_df.to_csv(output_path, sep='\t', index=False, columns=columns_to_output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--target_id_counts', required=True, help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')

    args = parser.parse_args()

    # Read the target_id_counts file
    target_id_counts_df = pd.read_csv(args.target_id_counts, sep='\t')

    distill_summary(args.combined_annotations, target_id_counts_df, args.output)
