import pandas as pd
import os
import argparse

def distill_summary(combined_annotations_path, distill_sheets_file, output_path):
    # Print the contents of distill_sheets_file
    print("distill_sheets_file:", distill_sheets_file)

    # Read the combined_annotations file
    combined_annotations_df = pd.read_csv(combined_annotations_path, sep='\t')

    # Identify columns ending in "_id" (excluding 'query_id')
    id_columns = [col for col in combined_annotations_df.columns if col.endswith("_id") and col != 'query_id']

    # Initialize an empty DataFrame to store the distill summary
    distill_summary_df = pd.DataFrame()

    # Read the distill_sheets_file and split paths
    with open(distill_sheets_file, 'r') as file:
        distill_sheets = file.read().split(',')

    # Process each distill sheet
    for distill_sheet in distill_sheets:
        distill_sheet = distill_sheet.strip()  # Remove leading/trailing whitespaces
        if not distill_sheet:
            continue  # Skip empty paths

        # Read the distill sheet
        distill_df = pd.read_csv(distill_sheet, sep='\t')

        # Identify the column to use as the sheet name
        sheet_name_column = [col for col in distill_df.columns if col.endswith("_ecosystem")]

        if not sheet_name_column:
            print(f"Error: No column ending with '_ecosystem' found in {distill_sheet}")
            continue

        sheet_name_column = sheet_name_column[0]

        # Use the values in the identified column as the sheet names
        sheet_names = distill_df[sheet_name_column].tolist()

        # Merge the distill sheet with the combined_annotations on 'gene_id' and identified "_id" columns
        merged_df = pd.merge(combined_annotations_df, distill_df, left_on=id_columns, right_on='gene_id', how='inner')

        # Append the merged DataFrame to the distill summary DataFrame
        distill_summary_df = pd.concat([distill_summary_df, merged_df])

    # Save the distill summary to the specified output path
    distill_summary_df.to_csv(output_path, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate genome summary from distill sheets and combined annotations.')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined_annotations.tsv file.')
    parser.add_argument('--distill_sheets_file', required=True, help='Path to the distill_sheets_file containing paths to TSV files.')
    parser.add_argument('--output', required=True, help='Path to the output genome_summary.tsv file.')

    args = parser.parse_args()

    distill_summary(args.combined_annotations, args.distill_sheets_file, args.output)
