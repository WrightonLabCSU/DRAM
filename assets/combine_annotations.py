import argparse
import pandas as pd
import logging
import os

# Configure the logger
logging.basicConfig(filename="logs/combine_annotations.log", level=logging.INFO, format='%(levelname)s: %(message)s')

def combine_annotations(annotation_files, output_file):
    # Create an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()

    # Process the input annotation files
    for i in range(0, len(annotation_files), 2):
        sample = annotation_files[i]
        file_path = annotation_files[i + 1]

        # Extract sample names from the input
        sample = sample.strip('\'"')

        # Check if the file exists
        if not os.path.exists(file_path):
            raise ValueError(f"Could not find file: {file_path}")

        # Read each annotation file
        logging.info(f"Processing annotation file: {file_path}")

        try:
            # Read CSV with custom separator and without specifying names
            annotation_data = pd.read_csv(file_path, sep=',', header=0)
        except FileNotFoundError:
            raise ValueError(f"Could not find file: {file_path}")

        # Print column names for each annotation file
        logging.info(f"Column names: {annotation_data.columns}")

        # Make the script case-insensitive when checking for 'query_id'
        query_id_col = [col for col in annotation_data.columns if col.lower() == 'query_id']
        if not query_id_col:
            raise ValueError(f"Column 'query_id' not found in the annotation file: {file_path}")
        query_id_col = query_id_col[0]

        # Merge dataframes based on 'query_id' column
        if combined_data.empty:
            combined_data = annotation_data.copy()
        else:
            common_cols = set(combined_data.columns) & set(annotation_data.columns)
            merge_cols = ['query_id'] + list(common_cols - {'query_id'})
            combined_data = pd.merge(combined_data, annotation_data, how='outer', on=merge_cols, suffixes=('_dbcan', '_kofam'))


        logging.info(f"Processed annotation file: {file_path} for sample: {sample}")

    # Rearrange the columns to match the desired order
    output_columns_order = ['query_id', 'sample'] + sorted([col for col in combined_data.columns if col not in ['query_id', 'sample']])
    combined_data = combined_data[output_columns_order]

    # Convert 'sample' column values to strings and join them with a semicolon
    sample_cols = [col for col in combined_data.columns if col.startswith('sample')]
    combined_data['sample'] = combined_data[sample_cols].apply(lambda row: "; ".join(map(str, filter(lambda x: pd.notna(x), row))), axis=1)

    # Save the combined data to the output file
    combined_data.to_csv(output_file, sep='\t', index=False)
    logging.info("Combining annotations completed.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine multiple annotation files into one.')
    parser.add_argument('--annotations', nargs='+', help='List of input annotation files', required=True)
    parser.add_argument('--output', help='Output file name', required=True)

    args = parser.parse_args()

    # Remove square brackets and extra commas
    args.annotations = [arg.strip("[],") for arg in args.annotations]

    combine_annotations(args.annotations, args.output)
