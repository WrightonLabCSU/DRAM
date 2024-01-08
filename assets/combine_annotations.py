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
            logging.warning(f"Could not find file: {file_path}")
            continue

        # Read each annotation file
        logging.info(f"Processing annotation file: {file_path}")

        try:
            # Read CSV with custom separator and without specifying names
            annotation_data = pd.read_csv(file_path, sep=',', header=0)
        except FileNotFoundError:
            logging.warning(f"Could not find file: {file_path}")
            continue

        # Make the script case-insensitive when checking for 'query_id'
        query_id_col = [col for col in annotation_data.columns if col.lower() == 'query_id']
        if not query_id_col:
            logging.warning(f"Column 'query_id' not found in the annotation file: {file_path}")
            continue
        query_id_col = query_id_col[0]

        # Merge the current annotation_data with combined_data based on 'query_id'
        combined_data = pd.concat([combined_data, annotation_data.set_index(query_id_col)], axis=1, join='outer')

        # Remove duplicate samples within the same row and separate them with a semicolon
        combined_data['sample'] = combined_data['sample'].apply(lambda x: "; ".join(list(set(x))))

        # Identify the column with the highest bitScore for each 'query_id'
        bitScore_cols = [col for col in combined_data.columns if col.endswith('_bitScore')]
        combined_data['best_bitScore_col'] = combined_data[bitScore_cols].idxmax(axis=1)

        # Extract the corresponding column information for the best bitScore
        combined_data['best_bitScore'] = combined_data.apply(lambda row: row[row['best_bitScore_col']], axis=1)

        # Rearrange the columns to match the desired order
        output_columns_order = [query_id_col, 'sample'] + list(combined_data.columns.difference([query_id_col, 'sample', 'best_bitScore_col', 'best_bitScore']))
        combined_data = combined_data[output_columns_order].reset_index()

        logging.info(f"Processed annotation file: {file_path}")

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
