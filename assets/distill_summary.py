import pandas as pd
import os
import argparse

def distill_summary(combined_annotations_path, distill_sheets_file, output_path):
    # Read the combined_annotations file
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')

    # Initialize an empty DataFrame to store the distill summary
    distill_summary_df = pd.DataFrame(columns=combined_annotations_df.columns)

    # Read the paths to distill sheets from the text file
    with open(distill_sheets_file, 'r') as file:
        distill_sheets = [line.strip() for line in file.readlines() if line.strip()]

    # Process each distill sheet
    for distill_sheet in distill_sheets:
        # Split the comma-separated paths into individual paths
        distill_sheet_paths = distill_sheet.split(',')

        for sheet_path in distill_sheet_paths:
            # Read the distill sheet
            distill_df = pd.read_csv(sheet_path, sep='\t')

            # Process each gene_id
            for gene_id in distill_df['gene_id']:
                # Filter rows with the current gene_id
                gene_rows = distill_df[distill_df['gene_id'] == gene_id]

                # Concatenate the rows to the distill summary DataFrame
                distill_summary_df = pd.concat([distill_summary_df, gene_rows], ignore_index=True)

    # Save the distill summary to the specified output path
    distill_summary_df.to_csv(output_path, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--distill_sheets', required=True, help='Path to the text file containing distill sheets.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')

    args = parser.parse_args()

    distill_summary(args.combined_annotations, args.distill_sheets, args.output)
