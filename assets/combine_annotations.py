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
    for annotation_file in annotation_files:
        # Extract sample name from the input file path
        sample = os.path.splitext(os.path.basename(annotation_file))[0]

        # Check if the file exists
        if not os.path.exists(annotation_file):
            raise ValueError(f"Could not find file: {annotation_file}")

        # Read each annotation file
        logging.info(f"Processing annotation file: {annotation_file}")

        try:
            # Read CSV without specifying names
            annotation_data = pd.read_csv(annotation_file, header=0)
        except FileNotFoundError:
            raise ValueError(f"Could not find file: {annotation_file}")

        # Print column names for each annotation file
        logging.info(f"Column names: {annotation_data.columns}")

        # Find the '*_bitScore' column dynamically
        bit_score_col = [col for col in annotation_data.columns if col.endswith('_bitScore')]
        if not bit_score_col:
            raise ValueError(f"No '*_bitScore' column found in the annotation file: {annotation_file}")
        bit_score_col = bit_score_col[0]

        # Find the row with the highest bitScore for each query_id
        annotation_data = annotation_data.sort_values(by=bit_score_col, ascending=False).drop_duplicates('query_id')

        # Merge dataframes based on 'query_id' column
        if combined_data.empty:
            combined_data = annotation_data.copy()
        else:
            combined_data = pd.concat([combined_data, annotation_data], ignore_index=True)

        logging.info(f"Processed annotation file: {annotation_file} for sample: {sample}")

    # Rearrange the columns to match the desired order
    output_columns_order = ['query_id'] + sorted([col for col in combined_data.columns if col != 'query_id'])
    combined_data = combined_data[output_columns_order]

    # Save the combined data to the output file
    combined_data.to_csv(output_file, sep='\t', index=False)
    logging.info("Combining annotations completed.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine multiple annotation files into one.')
    parser.add_argument('--annotations', nargs='+', help='List of input annotation files', required=True)
    parser.add_argument('--output', help='Output file name', required=True)

    args = parser.parse_args()

    combine_annotations(args.annotations, args.output)
