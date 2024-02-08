import pandas as pd
import argparse
import glob
import logging

logging.basicConfig(level=logging.DEBUG)


def is_null_content(file_path):
    with open(file_path, 'r') as file:
        content = file.read().strip()
    return content == "NULL"


def partial_match(gene_id, combined_column):
    return gene_id in combined_column


def distill_summary(combined_annotations_path, target_id_counts_df, output_path):
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')
    potential_gene_id_columns = [col for col in combined_annotations_df.columns if
                                 col.endswith('_id') and col != "query_id"]
    potential_ec_columns = [col for col in combined_annotations_df.columns if col.endswith('_EC')]
    distill_sheets = glob.glob('*_distill_sheet.tsv')
    distill_summary_df = pd.DataFrame()
    additional_columns = set()

    for distill_sheet in distill_sheets:
        if is_null_content(distill_sheet):
            logging.info(f"Skipping distill sheet '{distill_sheet}' as it contains 'NULL' content.")
            continue
        logging.info(f"Processing distill sheet: {distill_sheet}")
        distill_df = pd.read_csv(distill_sheet, sep='\t')

        # Collect additional columns that are not in the combined_annotations_df
        additional_columns.update(set(distill_df.columns) - set(combined_annotations_df.columns) - {'gene_id'})

        for common_gene_id_column in potential_gene_id_columns + potential_ec_columns:
            # Filter combined_annotations based on partial matching
            partial_match_indices = distill_df['gene_id'].apply(lambda x: partial_match(x, combined_annotations_df[common_gene_id_column]))
            partial_matched_combined_annotations = combined_annotations_df[partial_match_indices]

            # Merge the distill sheet with the filtered combined_annotations
            merged_df = pd.merge(
                partial_matched_combined_annotations,
                distill_df,
                left_on=[common_gene_id_column],
                right_on=['gene_id'],
                how='inner'
            )
            distill_summary_df = pd.concat([distill_summary_df, merged_df], ignore_index=True)

    distill_summary_df = pd.merge(distill_summary_df, target_id_counts_df, left_on=['gene_id'], right_on=['target_id'],
                                  how='left')
    # Remove the 'target_id' column to avoid duplicates
    distill_summary_df.drop('target_id', axis=1, inplace=True, errors='ignore')

    columns_to_output = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']
    columns_to_output.extend(additional_columns)
    # Make sure to remove any duplicates from the columns_to_output
    columns_to_output = list(dict.fromkeys(columns_to_output))

    bin_columns = [col for col in target_id_counts_df.columns if col.startswith('bin-')]
    columns_to_output.extend(bin_columns)

    deduplicated_df = distill_summary_df.drop_duplicates(subset=columns_to_output[:6]).copy()
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
