import pandas as pd
import argparse
import glob
import logging
import re
from pandas import concat

# Setup logging
logging.basicConfig(level=logging.DEBUG)

def is_null_content(file_path):
    """Check if the content of a file is NULL."""
    with open(file_path, 'r') as file:
        content = file.read().strip()
    return content == "NULL"

def is_partial_match(ec_number, partial_ec):
    """
    Check if the EC number starts with the given partial EC number and optionally followed by more subdivisions
    or a dash indicating unspecified subdivisions.

    Args:
        ec_number (str): The EC number to check.
        partial_ec (str): The partial EC number to match the start against.

    Returns:
        bool: True if the EC number starts with the partial EC number and optionally followed by more subdivisions or a dash, False otherwise.
    """
    if not isinstance(ec_number, str):
        return False

    # Build a regex pattern that starts with the partial_ec followed by any number of additional subdivisions
    # or a dash, which may be at the end or followed by further subdivisions
    pattern = re.compile(rf'^{re.escape(partial_ec)}(\.\d+)*(\.-)?$')
    return bool(pattern.match(ec_number))
def distill_summary(combined_annotations_path, target_id_counts_df, output_path):
    """
    Generate a genome summary from distill sheets and combined annotations.

    Args:
        combined_annotations_path (str): Path to the combined_annotations.tsv file.
        target_id_counts_df (pandas.DataFrame): DataFrame containing target ID counts.
        output_path (str): Path to the output genome_summary.tsv file.
    """
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')
    potential_gene_id_columns = [col for col in combined_annotations_df.columns if col.endswith('_id') and col != "query_id"]
    distill_sheets = glob.glob('*_distill_sheet.tsv')
    distill_summary_df = pd.DataFrame()
    has_target_id_column = False

    required_columns = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']
    # Initialize an empty list to store additional columns
    additional_columns = []

    for distill_sheet in distill_sheets:
        if is_null_content(distill_sheet):
            logging.info(f"Skipping distill sheet '{distill_sheet}' as it contains 'NULL' content.")
            continue
        logging.info(f"Processing distill sheet: {distill_sheet}")
        distill_df = pd.read_csv(distill_sheet, sep='\t')
        logging.debug(f"DataFrame from {distill_sheet}:")
        logging.debug(distill_df.head())

        for gene_id, row in distill_df.groupby('gene_id'):
            gene_description = row['gene_description'].iloc[0]  # Get the first value of 'gene_description' (assuming it's the same for all rows with the same 'gene_id')
            pathway = row['pathway'].iloc[0] if 'pathway' in row else None
            topic_ecosystem = row['topic_ecosystem'].iloc[0] if 'topic_ecosystem' in row else None
            category = row['category'].iloc[0] if 'category' in row else None
            subcategory = row['subcategory'].iloc[0] if 'subcategory' in row else None

            gene_id_found = False  # Flag to check if gene_id is found in any potential column or potential EC column

            # Inside the loop where potential gene ID columns are checked
            for col in potential_gene_id_columns:
                matched_indices = combined_annotations_df[col].str.contains('^' + re.escape(gene_id) + '$', na=False)
                if matched_indices.any():
                    gene_id_found = True
                    logging.debug(f"gene_id {gene_id} matched in column {col} with values:")
                    logging.debug(combined_annotations_df.loc[matched_indices, col].tolist())
                    
                    for combined_id in combined_annotations_df.loc[matched_indices, col]:
                        row_data = {
                            'gene_id': combined_id,
                            'gene_description': gene_description,
                            'pathway': pathway,
                            'category': category,
                            'subcategory': subcategory
                        }
                        # Include additional columns from the distill sheet
                        for additional_col in set(distill_df.columns) - set(combined_annotations_df.columns) - {'gene_id'}:
                            if additional_col in required_columns:
                                new_col_name = additional_col
                            else:
                                new_col_name = f"{additional_col}-{topic_ecosystem}"
                                additional_columns.append(new_col_name)
                            row_data[new_col_name] = row[additional_col].iloc[0] if additional_col in row else None
                        distill_summary_df = concat([distill_summary_df, pd.DataFrame([row_data])], ignore_index=True)
                    break


            # If gene_id is not found in any potential gene ID column, check potential EC columns
            if not gene_id_found:
                for col in combined_annotations_df.columns:
                    if col.endswith('_EC'):
                        for idx, ec_value in combined_annotations_df[col].iteritems():
                            if is_partial_match(ec_value, gene_id):
                                gene_id_found = True
                                logging.debug(f"Partial EC match found for gene_id {gene_id} in column {col}: {ec_value}")
                                
                                # Find associated gene_id values in the same row
                                associated_gene_ids = combined_annotations_df.loc[idx, [col.replace('_EC', '_id') for col in combined_annotations_df.columns if col.endswith('_id') and col != "query_id"]].tolist()
                                
                                for associated_id in associated_gene_ids:
                                    row_data = {
                                        'gene_id': associated_id,
                                        'gene_description': gene_description,
                                        'pathway': pathway,
                                        'category': category,
                                        'subcategory': subcategory,
                                        'associated_EC': gene_id  # Setting the associated_EC to the original gene_id
                                    }
                                    # Include additional columns from the distill sheet
                                    for additional_col in set(distill_df.columns) - set(combined_annotations_df.columns) - {'gene_id'}:
                                        if additional_col in required_columns:
                                            new_col_name = additional_col
                                        else:
                                            new_col_name = f"{additional_col}-{topic_ecosystem}"
                                            additional_columns.append(new_col_name)
                                        row_data[new_col_name] = row[additional_col].iloc[0] if additional_col in row else None
                                    distill_summary_df = concat([distill_summary_df, pd.DataFrame([row_data])], ignore_index=True)
                                break

    # Inside the distill_summary function after the for loop for processing distill sheets

    logging.debug("Distill summary DataFrame before merging:")
    logging.debug(distill_summary_df.head())

    # Check if 'gene_id' column exists in distill_summary_df
    if 'gene_id' not in distill_summary_df.columns:
        logging.error("No 'gene_id' column found in distill_summary_df. Aborting merge.")
        return

    # Inside the distill_summary function, before merging
    logging.debug("Distill summary DataFrame before merging:")
    logging.debug(distill_summary_df.head())

    # Check if 'target_id' column exists in target_id_counts_df
    if 'target_id' not in target_id_counts_df.columns:
        logging.error("No 'target_id' column found in target_id_counts_df. Aborting merge.")
        return

    # Merge with target_id_counts_df using 'target_id' column
    distill_summary_df = pd.merge(distill_summary_df, target_id_counts_df, left_on=['gene_id'], right_on=['target_id'], how='left')

    logging.debug("Distill summary DataFrame after merging:")
    logging.debug(distill_summary_df.head())

    # Drop the 'target_id' column if it was added from target_id_counts_df
    if 'target_id' in distill_summary_df.columns:
        distill_summary_df.drop('target_id', axis=1, inplace=True)

    # Drop rows with null gene_id
    deduplicated_df = distill_summary_df[~distill_summary_df['gene_id'].isnull()]

    logging.debug("Deduplicated DataFrame:")
    logging.debug(deduplicated_df.head())

    # Inside the distill_summary function after processing distill sheets

    # Add this line after reading distill_df
    logging.debug(f"DataFrame from {distill_sheet}:")
    logging.debug(distill_df.head())

    # Add this line inside the loop where potential gene ID columns are checked
    logging.debug(f"Checking potential gene ID columns for gene_id: {gene_id}")

    # Add this line after the loop where potential gene ID columns are checked
    logging.debug("Summary DataFrame after processing distill sheet:")
    logging.debug(distill_summary_df.head())

    # Write the resulting DataFrame to the output file
    deduplicated_df.to_csv(output_path, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--target_id_counts', required=True, help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')
    args = parser.parse_args()

    # Load target_id_counts DataFrame
    target_id_counts_df = pd.read_csv(args.target_id_counts, sep='\t')

    # Call distill_summary function with provided arguments
    distill_summary(args.combined_annotations, target_id_counts_df, args.output)