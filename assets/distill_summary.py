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

    for distill_sheet in distill_sheets:
        if is_null_content(distill_sheet):
            logging.info(f"Skipping distill sheet '{distill_sheet}' as it contains 'NULL' content.")
            continue
        logging.info(f"Processing distill sheet: {distill_sheet}")
        distill_df = pd.read_csv(distill_sheet, sep='\t')

        additional_cols = [col for col in distill_df.columns if col not in ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']]

        for gene_id, row in distill_df.groupby('gene_id'):
            gene_description = row['gene_description'].iloc[0]  # Get the first value of 'gene_description' (assuming it's the same for all rows with the same 'gene_id')
            pathway = row.get('pathway', None)
            topic_ecosystem = row.get('topic_ecosystem', None)
            category = row.get('category', None)
            subcategory = row.get('subcategory', None)

            gene_id_found = False

            for col in potential_gene_id_columns:
                matched_indices = combined_annotations_df[col].str.contains('^' + re.escape(gene_id) + '$', na=False)
                if matched_indices.any():
                    gene_id_found = True
                    
                    for combined_id in combined_annotations_df.loc[matched_indices, col]:
                        # Create a dictionary with row data including additional columns
                        row_data = {
                            'gene_id': combined_id,
                            'gene_description': gene_description,
                            'pathway': pathway,
                            'topic_ecosystem': topic_ecosystem,
                            'category': category,
                            'subcategory': subcategory
                        }
                        for additional_col in additional_cols:
                            row_data[additional_col] = row.get(additional_col, None)
                        if 'target_id' in row:
                            has_target_id_column = True
                        distill_summary_df = concat([distill_summary_df, pd.DataFrame([row_data])], ignore_index=True)
                    break

            if not gene_id_found:
                for col in combined_annotations_df.columns:
                    if col.endswith('_EC'):
                        for idx, ec_value in combined_annotations_df[col].items():  # Change iteritems() to items()
                            if is_partial_match(ec_value, gene_id):
                                gene_id_found = True
                                break

            if not gene_id_found:
                continue

            for col in combined_annotations_df.columns:
                if col.endswith('_EC'):
                    for idx, ec_value in combined_annotations_df[col].items():  # Change iteritems() to items()
                        ec_segments = str(ec_value).split(';')
                        for segment in ec_segments:
                            segment = segment.strip()
                            if is_partial_match(segment, gene_id):
                                associated_ec = segment

                                row_data = {
                                    'gene_id': None,
                                    'gene_description': gene_description,
                                    'pathway': pathway,
                                    'topic_ecosystem': topic_ecosystem,
                                    'category': category,
                                    'subcategory': subcategory,
                                    'associated_EC': associated_ec
                                }
                                for id_col in combined_annotations_df.filter(like='_id').columns:
                                    if id_col != 'query_id':
                                        gene_id_value = combined_annotations_df.at[idx, id_col]
                                        if gene_id_value:
                                            row_data['gene_id'] = gene_id_value
                                            distill_summary_df = pd.concat([distill_summary_df, pd.DataFrame([row_data])], ignore_index=True)
                                break

    distill_summary_df = pd.merge(distill_summary_df, target_id_counts_df, left_on=['gene_id'], right_on=['target_id'], how='left')
    
    if not has_target_id_column:
        distill_summary_df.drop('target_id', axis=1, inplace=True, errors='ignore')

    required_columns = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']
    if 'associated_EC' in distill_summary_df.columns:
        required_columns.append('associated_EC')
    additional_columns = [col for col in distill_summary_df.columns if col not in required_columns and col not in target_id_counts_df.columns]
    columns_to_output = required_columns + list(set(additional_columns) - {'associated_EC'}) + list(target_id_counts_df.columns)

    for col in columns_to_output:
        if col not in distill_summary_df.columns:
            distill_summary_df[col] = None

    try:
        # Deduplicate based on index
        deduplicated_df = distill_summary_df[~distill_summary_df.duplicated(keep='first')].copy()

        deduplicated_df = deduplicated_df[~deduplicated_df['gene_id'].isnull()]

        deduplicated_df.to_csv(output_path, sep='\t', index=False)
    except TypeError as e:
        print("TypeError occurred. Unhashable values:")
        for col, vals in distill_summary_df.iteritems():
            try:
                pd.Series(vals).unique()
            except TypeError:
                print(f"Column '{col}' contains unhashable values.")
        raise e

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--target_id_counts', required=True, help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')
    args = parser.parse_args()
    
    target_id_counts_df = pd.read_csv(args.target_id_counts, sep='\t')
    distill_summary(args.combined_annotations, target_id_counts_df, args.output)
