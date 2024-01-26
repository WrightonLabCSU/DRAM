import pandas as pd
import os
import argparse

def distill_summary(combined_annotations_path, distill_sheets, output_path):
    # Print the contents of distill_sheets
    print("distill_sheets:", distill_sheets)

    # Read the combined_annotations file
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')

    # Initialize an empty DataFrame to store the distill summary
    distill_summary_df = pd.DataFrame(columns=combined_annotations_df.columns)

    # Process each distill sheet
    for distill_sheet in distill_sheets:
        # Read the distill sheet
        distill_df = pd.read_csv(distill_sheet, sep='\t')

        # Use the values in the 'topic_ecosystem' column as the sheet names
        topic_ecosystem = distill_df['topic_ecosystem'].iloc[0]

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

    # Read the paths to distill sheets from the text file
    with open(args.distill_sheets, 'r') as file:
        distill_sheets = [line.strip() for line in file.readlines() if line.strip()]

    distill_summary(args.combined_annotations, distill_sheets, args.output)
