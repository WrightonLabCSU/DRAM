import argparse
import pandas as pd

def distill_summary(combined_annotations_path, distill_sheets, output_path):
    # Read input files
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t')

    # Check if the required columns are present in the combined_annotations
    required_columns = ['query_id', 'sample']
    missing_columns = set(required_columns) - set(combined_annotations.columns)

    if missing_columns:
        raise ValueError(f"The following required columns are missing in combined_annotations: {', '.join(missing_columns)}")

    # Initialize the output DataFrame with query_id, sample, and gene_id columns
    distill_summary_df = pd.DataFrame(columns=['query_id', 'sample', 'gene_id'])

    # Iterate through distill_sheets
    for distill_sheet in distill_sheets:
        # Read the distill_sheet
        distill_data = pd.read_csv(distill_sheet, sep='\t')

        # Check if the required columns are present in the distill_data
        required_columns_distill = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']
        missing_columns_distill = set(required_columns_distill) - set(distill_data.columns)

        if missing_columns_distill:
            raise ValueError(f"The following required columns are missing in {os.path.basename(distill_sheet)}: {', '.join(missing_columns_distill)}")

        # Merge combined_annotations with distill_data on 'gene_id'
        merged_data = pd.merge(combined_annotations, distill_data, on='gene_id', how='inner')

        # Add the relevant information to the output DataFrame
        distill_summary_df = pd.concat([
            distill_summary_df,
            merged_data[['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']]
        ])

    # Write the output to a TSV file
    distill_summary_df.to_csv(output_path, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--distill_sheets', nargs='+', required=True, help='List of paths to distill sheets')
    parser.add_argument('--output', required=True, help='Path to the output file')

    args = parser.parse_args()

    # Call the distill_summary function with provided arguments
    distill_summary(args.combined_annotations, args.distill_sheets, args.output)
