import pandas as pd
import os
import argparse

def distill_summary(combined_annotations_path, distill_sheets_file, output_path):
    # Read the combined_annotations file
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')

    # Parse the distill_sheets_file and process each distill sheet
    with open(distill_sheets_file, 'r') as file:
        distill_sheets = file.read().strip().split(',')
    
    # Initialize an empty DataFrame to store the distill summary
    distill_summary_df = pd.DataFrame()

    # Process the first distill sheet (for testing purposes)
    distill_sheet = distill_sheets[0]
    distill_df = pd.read_csv(distill_sheet, sep='\t')

    # Process each distill sheet
    for distill_sheet in distill_sheets:
        # Read the distill sheet
        distill_df = pd.read_csv(distill_sheet, sep='\t')

        # Extract the column name with "_id" (excluding query_id)
        id_column = [col for col in combined_annotations_df.columns if col.endswith('_id') and col != 'query_id'][0]

        # Merge the distill sheet with the combined_annotations on the appropriate ID column
        merged_df = pd.merge(combined_annotations_df, distill_df, left_on=id_column, right_on='gene_id', how='inner')

        # Append the merged DataFrame to the distill summary DataFrame
        distill_summary_df = pd.concat([distill_summary_df, merged_df])

    # Save the distill summary to the specified output path with only the required columns
    required_columns = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']
    distill_summary_df[required_columns].to_csv(output_path, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--distill_sheets_file', required=True, help='Path to the distill_sheets file containing paths to TSVs.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')

    args = parser.parse_args()

    distill_summary(args.combined_annotations, args.distill_sheets_file, args.output)
