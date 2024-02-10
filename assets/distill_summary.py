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

        for gene_id, row in distill_df.groupby('gene_id'):
            gene_description = row['gene_description'].iloc[0]  # Get the first value of 'gene_description' (assuming it's the same for all rows with the same 'gene_id')
            pathway = row['pathway'].iloc[0] if 'pathway' in row else None
            topic_ecosystem = row['topic_ecosystem'].iloc[0] if 'topic_ecosystem' in row else None
            category = row['category'].iloc[0] if 'category' in row else None
            subcategory = row['subcategory'].iloc[0] if 'subcategory' in row else None

            gene_id_found = False  # Flag to check if gene_id is found in any potential column or potential EC column

            # Check potential gene ID columns
            for col in potential_gene_id_columns:
                matched_indices = combined_annotations_df[col].str.contains('^' + re.escape(gene_id) + '$', na=False)
                if matched_indices.any():
                    gene_id_found = True
                    print(f"gene_id {gene_id} matched in column {col} with values:")
                    print(combined_annotations_df.loc[matched_indices, col].tolist())
                    
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
                            new_col_name = f"{additional_col}-{topic_ecosystem}"
                            if new_col_name not in required_columns:
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
                                print(f"Partial EC match found for gene_id {gene_id} in column {col}: {ec_value}")
                                
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
                                        new_col_name = f"{additional_col}-{topic_ecosystem}"
                                        if new_col_name not in required_columns:
                                            additional_columns.append(new_col_name)
                                            row_data[new_col_name] = row[additional_col].iloc[0] if additional_col in row else None
                                    distill_summary_df = concat([distill_summary_df, pd.DataFrame([row_data])], ignore_index=True)
                                break



    distill_summary_df = pd.merge(distill_summary_df, target_id_counts_df, left_on=['gene_id'], right_on=['target_id'], how='left')
    
    if not has_target_id_column:
        distill_summary_df.drop('target_id', axis=1, inplace=True)

    # Drop rows with null gene_id
    deduplicated_df = distill_summary_df[~distill_summary_df['gene_id'].isnull()]

    deduplicated_df.to_csv(output_path, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--target_id_counts', required=True, help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')
    args = parser.parse_args()

    target_id_counts_df = pd.read_csv(args.target_id_counts, sep='\t')
    distill_summary(args.combined_annotations, target_id_counts_df, args.output)
