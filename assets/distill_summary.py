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

def is_partial_match(gene_id, ec_number):
    """
    Check if the gene ID partially matches the given EC number.

    Args:
        gene_id (str): The gene ID to check.
        ec_number (str): The EC number to match against.

    Returns:
        bool: True if there is a partial match, False otherwise.
    """
    if not isinstance(gene_id, str) or not isinstance(ec_number, str):
        return False
    
    gene_id_segments = gene_id.split('.')
    ec_number_segments = ec_number.split('.')
    
    for i in range(len(ec_number_segments)):
        if i >= len(gene_id_segments) or gene_id_segments[i] != ec_number_segments[i]:
            return False
    
    return True

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

        for gene_id, row in distill_df.groupby('gene_id'):
            gene_description = row['gene_description'].iloc[0]  # Get the first value of 'gene_description' (assuming it's the same for all rows with the same 'gene_id')
            pathway = row['pathway'].iloc[0] if 'pathway' in row else None
            topic_ecosystem = row['topic_ecosystem'].iloc[0] if 'topic_ecosystem' in row else None
            category = row['category'].iloc[0] if 'category' in row else None
            subcategory = row['subcategory'].iloc[0] if 'subcategory' in row else None
            level = row['level'].iloc[0] if 'level' in row else None

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
                            'topic_ecosystem': topic_ecosystem,
                            'category': category,
                            'subcategory': subcategory,
                            'level': level
                        }
                        for additional_col in set(distill_df.columns) - set(combined_annotations_df.columns) - {'gene_id'}:
                            if additional_col == 'target_id':
                                has_target_id_column = True
                            row_data[additional_col] = row[additional_col].iloc[0] if additional_col in row else None
                        distill_summary_df = concat([distill_summary_df, pd.DataFrame([row_data])], ignore_index=True)
                    break

            # If gene_id is not found in any potential gene ID column, check partial matches with EC numbers
            if not gene_id_found:
                for ec_col in combined_annotations_df.columns:
                    if ec_col.endswith('_EC'):
                        for ec_value in combined_annotations_df[ec_col]:
                            if is_partial_match(gene_id, ec_value):
                                gene_id_found = True
                                print(f"Partial EC match found for gene_id {gene_id} in column {ec_col}: {ec_value}")
                                break

            # If gene_id is still not found, skip processing this gene_id
            if not gene_id_found:
                continue

            # Process associated EC values
            for col in combined_annotations_df.columns:
                if col.endswith('_EC'):
                    for idx, ec_value in combined_annotations_df[col].iteritems():
                        if gene_id in str(ec_value):
                            ec_segments = str(ec_value).split(';')
                            for segment in ec_segments:
                                if gene_id in segment:
                                    associated_ec = segment.strip()

                                    row_data = {
                                        'gene_id': None,
                                        'gene_description': gene_description,
                                        'pathway': pathway,
                                        'topic_ecosystem': topic_ecosystem,
                                        'category': category,
                                        'subcategory': subcategory,
                                        'level': level
                                    }
                                    if associated_ec:
                                        row_data['associated_EC'] = associated_ec
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

    required_columns = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory', 'level']
    if 'associated_EC' in distill_summary_df.columns:
        required_columns.append('associated_EC')
    additional_columns = [col for col in distill_summary_df.columns if col not in required_columns and col not in target_id_counts_df.columns]
    columns_to_output = required_columns + list(set(additional_columns) - {'associated_EC'}) + list(target_id_counts_df.columns)

    for col in columns_to_output:
        if col not in distill_summary_df.columns:
            distill_summary_df[col] = None

    deduplicated_df = distill_summary_df.drop_duplicates(subset=required_columns, ignore_index=True).copy()


    deduplicated_df.to_csv(output_path, sep='\t', index=False, columns=columns_to_output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--target_id_counts', required=True, help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')
    args = parser.parse_args()
    
    target_id_counts_df = pd.read_csv(args.target_id_counts, sep='\t')
    distill_summary(args.combined_annotations, target_id_counts_df, args.output)
