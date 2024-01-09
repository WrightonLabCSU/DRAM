import argparse
import pandas as pd
import logging

# Configure the logger
logging.basicConfig(filename="logs/combine_annotations.log", level=logging.INFO, format='%(levelname)s: %(message)s')

def combine_annotations(annotation_files, output_file):
    # Create an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()

    # Process the input annotation files
    for i in range(0, len(annotation_files), 2):
        # Extract sample name from the input annotation_files
        sample = annotation_files[i].split(',')[0]
        # Remove quotes from the sample name if present
        sample = sample.strip("'")

        file_path = annotation_files[i + 1]

        # Read each annotation file
        logging.info(f"Processing annotation file: {file_path}")

        try:
            annotation_data = pd.read_csv(file_path, sep=',', quotechar='"')
        except FileNotFoundError:
            raise ValueError(f"Could not find file: {file_path}")

        # Merge based on 'query_id' using an outer join
        combined_data = pd.merge(combined_data, annotation_data, on='query_id', how='outer')

        logging.info(f"Processed annotation file: {file_path}")

    # Reorder columns alphabetically
    combined_data = combined_data.reindex(sorted(combined_data.columns), axis=1)

    # Extract 'query_id' and 'sample' columns to come first
    col_order = ['query_id', 'sample'] + sorted([col for col in combined_data.columns if col not in ['query_id', 'sample']])

    # Reorder columns in the desired order
    combined_data = combined_data[col_order]

    # Save the combined data to the output file
    combined_data.to_csv(output_file, sep='\t', index=False)
    logging.info("Combining annotations completed.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine multiple annotation files into one.')
    parser.add_argument('--annotations', nargs='+', help='List of input annotation files and sample names', required=True)
    parser.add_argument('--output', help='Output file name', required=True)

    args = parser.parse_args()

    # Remove square brackets and extra commas
    args.annotations = [arg.strip("[],") for arg in args.annotations]

    combine_annotations(args.annotations, args.output)
