import argparse
import pandas as pd
import logging
import os

# Configure the logger
logging.basicConfig(filename="logs/combine_annotations.log", level=logging.INFO, format='%(levelname)s: %(message)s')

def combine_annotations(annotation_files, output_file):
    # Create an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()

    # Create a dictionary to store values for each unique query_id
    data_dict = {}

    # Set to store all unique column names
    all_columns = set()

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

        # Update set of all columns
        all_columns.update(annotation_data.columns)

        for index, row in annotation_data.iterrows():
            query_id = row[query_id_col]

            # Check if query_id already exists in the dictionary
            if query_id in data_dict:
                # Choose the row with the highest bitScore
                db_name = os.path.basename(file_path).split('_')[0]
                if f'{db_name}_bitScore' in row.index and f'{db_name}_bitScore' in data_dict[query_id].index:
                    if row[f'{db_name}_bitScore'] > data_dict[query_id][f'{db_name}_bitScore']:
                        data_dict[query_id] = row
                    # Append the sample to the list
                    if sample not in data_dict[query_id]['sample']:
                        data_dict[query_id]['sample'].append(sample)
            else:
                # Create a new entry in the dictionary
                data_dict[query_id] = row
                data_dict[query_id]['sample'] = [sample]

            logging.info(f"Processed query_id: {query_id} for sample: {sample}")

    # Create a DataFrame from the dictionary
    combined_data = pd.DataFrame.from_dict(data_dict, orient='index')

    # Rearrange the columns to match the desired order
    output_columns_order = [query_id_col, 'sample'] + list(all_columns.difference([query_id_col, 'sample']))
    combined_data = combined_data[output_columns_order]

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
