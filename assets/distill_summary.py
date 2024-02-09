import pandas as pd
import argparse
import glob
import logging
import re
from pandas import concat

# Setup logging
logging.basicConfig(level=logging.DEBUG)

def is_null_content(file_path):
    with open(file_path, 'r') as file:
        content = file.read().strip()
    return content == "NULL"

def is_partial_match(ec_number, partial_ec):
    # Check if ec_number is a string or bytes-like object
    if not isinstance(ec_number, (str, bytes)):
        return False
    
    # Escape the dots in the partial EC number to match them as literals
    partial_ec_escaped = re.escape(partial_ec)
    # Construct a regular expression pattern to match the partial EC number at the start
    # and ensure it is followed by a dot and more digits
    pattern = re.compile(rf'^{partial_ec_escaped}(\.\d+)+$')
    # Check if the EC number matches the pattern
    return bool(pattern.match(ec_number))


def distill_summary(combined_annotations_path, target_id_counts_df, output_path):
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')
    potential_gene_id_columns = [col for col in combined_annotations_df.columns if col.endswith('_id') and col != "query_id"]
    distill_sheets = glob.glob('*_distill_sheet.tsv')
    distill_summary_df = pd.DataFrame()
    additional_columns = set()
    has_target_id_column = False

    for distill_sheet in distill_sheets:
        if is_null_content(distill_sheet):
            logging.info(f"Skipping distill sheet '{distill_sheet}' as it contains 'NULL' content.")
            continue
        logging.info(f"Processing distill sheet: {distill_sheet}")
        distill_df = pd.read_csv(distill_sheet, sep='\t')

        for index, row in distill_df.iterrows():
            gene_id = row['gene_id']
            gene_description = row['gene_description']
            pathway = row.get('pathway', None)
            topic_ecosystem = row.get('topic_ecosystem', None)
            category = row.get('category', None)
            subcategory = row.get('subcategory', None)

            # First, try matching gene_id against potential_gene_id_columns
            for col in potential_gene_id_columns:
                matched_indices = combined_annotations_df[col].str.contains(gene_id, na=False)
                if matched_indices.any():
                    for combined_id in combined_annotations_df.loc[matched_indices, col]:
                        row_data = {
                            'gene_id': combined_id,
                            'gene_description': gene_description,
                            'pathway': pathway,
                            'topic_ecosystem': topic_ecosystem,
                            'category': category,
                            'subcategory': subcategory
                            # Note: 'associated_EC' not added here as it's specific to EC number associations
                        }
                        # Add additional columns from distill sheet
                        for additional_col in set(distill_df.columns) - set(combined_annotations_df.columns) - {'gene_id'}:
                            if additional_col == 'target_id':
                                has_target_id_column = True
                            row_data[additional_col] = row.get(additional_col, None)
                        distill_summary_df = concat([distill_summary_df, pd.DataFrame([row_data])], ignore_index=True)
                    break  # Break after matching to avoid processing the same gene_id against multiple columns
            else:
                for col in combined_annotations_df.columns:
                    if col.endswith('_EC'):
                        for index, ec_value in combined_annotations_df[col].iteritems():
                            # Check for partial matches within the cell and capture the full EC number
                            if gene_id in str(ec_value):
                                # Split the cell by semicolon and look for the gene_id within each segment
                                ec_segments = str(ec_value).split(';')
                                for segment in ec_segments:
                                    if gene_id in segment:
                                        associated_ec = segment.strip()  # This should be the full EC number
                                        print(f"Match found for gene_id {gene_id} in column {col}: {associated_ec}")

                                        # Now we can gather corresponding gene_ids from all _id columns in the same row
                                        row_data = {
                                            'gene_id': None,  # Placeholder for now, will be updated below
                                            'gene_description': gene_description,
                                            'pathway': pathway,
                                            'topic_ecosystem': topic_ecosystem,
                                            'category': category,
                                            'subcategory': subcategory,
                                            'associated_EC': associated_ec  # Include the full EC number
                                        }
                                        for id_col in combined_annotations_df.filter(like='_id').columns:
                                            if id_col != 'query_id':  # Skip query_id column
                                                gene_id_value = combined_annotations_df.at[index, id_col]
                                                if gene_id_value:  # Ensure the value is not empty
                                                    row_data['gene_id'] = gene_id_value
                                                    distill_summary_df = pd.concat([distill_summary_df, pd.DataFrame([row_data])], ignore_index=True)
                                        break  # Stop searching for further matches in this cell



    # Merge distill_summary_df with target_id_counts_df
    distill_summary_df = pd.merge(distill_summary_df, target_id_counts_df, left_on=['gene_id'], right_on=['target_id'],
                                  how='left')
    
    # Remove the 'target_id' column if it wasn't added as an additional column
    if not has_target_id_column:
        distill_summary_df.drop('target_id', axis=1, inplace=True, errors='ignore')

    # Define columns to output
    required_columns = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory', 'associated_EC']
    additional_columns = [col for col in distill_summary_df.columns if col not in required_columns and col not in target_id_counts_df.columns]

    columns_to_output = required_columns + list(set(additional_columns) - {'associated_EC'}) + list(target_id_counts_df.columns)

    # Ensure all required columns are present
    for col in columns_to_output:
        if col not in distill_summary_df.columns:
            distill_summary_df[col] = None

    # Drop duplicates based on subset of columns
    deduplicated_df = distill_summary_df.drop_duplicates(subset=required_columns, ignore_index=True).copy()

    # Write output to file
    deduplicated_df.to_csv(output_path, sep='\t', index=False, columns=columns_to_output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--target_id_counts', required=True, help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')

    args = parser.parse_args()
    target_id_counts_df = pd.read_csv(args.target_id_counts, sep='\t')
    distill_summary(args.combined_annotations, target_id_counts_df, args.output)

