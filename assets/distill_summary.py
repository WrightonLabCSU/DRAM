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
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')
    potential_gene_id_columns = [col for col in combined_annotations_df.columns if col.endswith('_id') and col != "query_id"]
    distill_sheets = glob.glob('*_distill_sheet.tsv')
    distill_summary_df = pd.DataFrame()

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

            # Check potential_gene_id_columns first
            for col in combined_annotations_df.columns:
                if col.endswith('_id') and col != "query_id":  # Exclude query_id
                    matched_indices = combined_annotations_df[col].str.contains(gene_id, na=False)
                    if matched_indices.any():
                        combined_ids = gene_id  # Use the matched gene_id directly
                        distill_summary_df = distill_summary_df.append({'gene_id': combined_ids,
                                                                        'gene_description': gene_description,
                                                                        'pathway': pathway,
                                                                        'topic_ecosystem': topic_ecosystem,
                                                                        'category': category,
                                                                        'subcategory': subcategory},
                                                                       ignore_index=True)
                        break
            else:
                # Check potential_ec_columns
                for col in combined_annotations_df.columns:
                    if col.endswith('_EC'):
                        matched_indices = combined_annotations_df[col].str.contains(gene_id, na=False)
                        if matched_indices.any():
                            combined_ids = '; '.join(combined_annotations_df.loc[matched_indices, col.replace('_EC', '_id')])
                            distill_summary_df = distill_summary_df.append({'gene_id': combined_ids,
                                                                            'gene_description': gene_description + '; ' + gene_id,
                                                                            'pathway': pathway,
                                                                            'topic_ecosystem': topic_ecosystem,
                                                                            'category': category,
                                                                            'subcategory': subcategory},
                                                                           ignore_index=True)
                            break
                else:
                    logging.info(f"No match found for gene ID {gene_id}.")

    # Merge distill_summary_df with target_id_counts_df
    distill_summary_df = pd.merge(distill_summary_df, target_id_counts_df, left_on=['gene_id'], right_on=['target_id'],
                                  how='left')
    
    # Remove the 'target_id' column to avoid duplicates
    distill_summary_df.drop('target_id', axis=1, inplace=True, errors='ignore')

    # Define columns to output
    columns_to_output = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']

    # Add bin columns from target_id_counts_df
    bin_columns = [col for col in target_id_counts_df.columns if col.startswith('bin-')]
    columns_to_output.extend(bin_columns)

    # Ensure all required columns are present
    for col in columns_to_output:
        if col not in distill_summary_df.columns:
            distill_summary_df[col] = None

    # Drop duplicates based on subset of columns
    deduplicated_df = distill_summary_df.drop_duplicates(subset=columns_to_output[:6], ignore_index=True).copy()

    # Write output to file
    deduplicated_df.to_csv(output_path, sep='\t', index=False, columns=columns_to_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True,
                        help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--target_id_counts', required=True,
                        help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--output', required=True,
                        help='Path to the output genome_summary.tsv file.')

    args = parser.parse_args()
    target_id_counts_df = pd.read_csv(args.target_id_counts, sep='\t')
    distill_summary(args.combined_annotations, target_id_counts_df, args.output)
