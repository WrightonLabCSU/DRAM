import argparse
import pandas as pd
import logging

# Configure the logger
logging.basicConfig(filename="logs/combine_annotations.log", level=logging.INFO, format='%(levelname)s: %(message)s')

def identify_columns(columns, suffix):
    matching_columns = [col for col in columns if col.lower().endswith(suffix.lower())]
    if not matching_columns:
        raise ValueError(f"No column ending with '{suffix}' found.")
    return matching_columns[0]

def combine_annotations(annotation_files, output_file):
    # Create an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()

    # Process the input annotation files
    for i in range(0, len(annotation_files), 2):
        sample = annotation_files[i]
        file_path = annotation_files[i + 1]

        # Extract sample names from the input
        sample = sample.split('\'')[1] if '\'' in sample else sample

        # Read each annotation file
        logging.info(f"Processing annotation file: {file_path}")

        try:
            annotation_data = pd.read_csv(file_path, sep=',', quotechar='"')
        except FileNotFoundError:
            raise ValueError(f"Could not find file: {file_path}")

        # Identify columns based on suffixes
        target_id_col = identify_columns(annotation_data.columns, "_id")
        score_rank_col = identify_columns(annotation_data.columns, "_score_rank")
        bitScore_col = identify_columns(annotation_data.columns, "_bitScore")

        # Merge additional columns with the existing columns (excluding 'query_id')
        additional_cols = annotation_data.columns.difference(['query_id'])
        cols_to_merge = ['query_id', 'sample'] + additional_cols.tolist()

        # Merge data based on 'query_id'
        combined_data = pd.merge(combined_data, annotation_data[cols_to_merge], on='query_id', how='outer')

        logging.info(f"Processed annotation file: {file_path}")

    # Reorder columns alphabetically
    combined_data = combined_data.reindex(sorted(combined_data.columns), axis=1)

    # Remove duplicate samples within the same row and separate them with a semicolon
    combined_data['sample'] = combined_data['sample'].apply(lambda x: "; ".join(list(set(x))))

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
