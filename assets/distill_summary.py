import pandas as pd
import os
import argparse

def distill_summary(combined_annotations_path, distill_sheets_file, output_path):
    # Read the combined_annotations file
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')

    # Read the distill_sheets_file and split paths
    with open(distill_sheets_file, 'r') as file:
        distill_sheets = file.read().strip().split(',')

    # Print the contents of distill_sheets
    print("distill_sheets:", distill_sheets)

    # Initialize an empty DataFrame to store the distill summary
    distill_summary_df = pd.DataFrame()

    # Process each distill sheet
    for distill_sheet in distill_sheets:
        # Read the distill sheet
        distill_df = pd.read_csv(distill_sheet.strip(), sep='\t')

        # Use the values in the 'topic_ecosystem' column as the sheet names
        topic_ecosystem = distill_df['topic_ecosystem'].iloc[0]

        # Merge the distill sheet with the combined_annotations on 'gene_id'
        merged_df = pd.merge(combined_annotations_df, distill_df, on='gene_id', how='inner')

        # Append the merged DataFrame to the distill summary DataFrame
        distill_summary_df = pd.concat([distill_summary_df, merged_df])

    # Save the distill summary to the specified output path
    distill_summary_df.to_csv(output_path, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--distill_sheets_file', required=True, help='Path to the text file containing paths to distill sheets.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')

    args = parser.parse_args()

    distill_summary(args.combined_annotations, args.distill_sheets_file, args.output)
