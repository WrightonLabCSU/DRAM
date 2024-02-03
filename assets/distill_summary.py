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
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')
    potential_gene_id_columns = [col for col in combined_annotations_df.columns if col.endswith('_id') and col != "query_id"]
    distill_sheets = glob.glob('*_distill_sheet.tsv')
    distill_summary_df = pd.DataFrame()
    sample_columns = target_id_counts_df.columns[target_id_counts_df.dtypes == 'int64']
    sample_names = sample_columns.tolist()

    print("Columns in target_id_counts:", sample_names)

    for distill_sheet in distill_sheets:
        if is_null_content(distill_sheet):
            logging.info(f"Skipping distill sheet '{distill_sheet}' as it contains 'NULL' content.")
            continue
        logging.info(f"Processing distill sheet: {distill_sheet}")
        distill_df = pd.read_csv(distill_sheet, sep='\t')
        logging.info(f"Column names of distill sheet: {distill_df.columns}")

        merged_data_for_current_gene_id = pd.DataFrame()

        for common_gene_id_column in potential_gene_id_columns:
            try:
                merged_df = pd.merge(
                    combined_annotations_df,
                    distill_df,
                    left_on=[common_gene_id_column],
                    right_on=['gene_id'],
                    how='inner'
                )
                merged_data_for_current_gene_id = pd.concat([merged_data_for_current_gene_id, merged_df])
            except KeyError:
                logging.warning(f"'{common_gene_id_column}' column not found in '{distill_sheet}'. Skipping merge for this column.")

        distill_summary_df = pd.concat([distill_summary_df, merged_data_for_current_gene_id])

    if "potential_amg" in distill_df.columns:
        distill_summary_df = pd.merge(
            distill_summary_df,
            target_id_counts_df,
            left_on=['gene_id'],
            right_on=['target_id'],
            how='left'
        )

        # Update sample columns with actual counts after merge
        for sample_name in sample_names:
            distill_summary_df[sample_name] = distill_summary_df[sample_name].fillna(0)

    deduplicated_df = distill_summary_df.drop_duplicates(subset=['gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory'])
    deduplicated_df.to_csv(output_path, sep='\t', index=False)

    # Debugging: Print final dataframe shape and sample columns
    print("Final DataFrame shape:", deduplicated_df.shape)
    print("Sample columns in final DataFrame:", deduplicated_df[sample_names].head())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--target_id_counts', required=True, help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')

    args = parser.parse_args()

    # Read the target_id_counts file
    target_id_counts_df = pd.read_csv(args.target_id_counts, sep='\t')

    distill_summary(args.combined_annotations, target_id_counts_df, args.output)
