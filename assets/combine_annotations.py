import argparse
import pandas as pd
import logging

# Configure the logger
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

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
        sample = sample.split('\'')[1] if '\'' in sample else sample

        # Read each annotation file
        logging.info(f"Processing annotation file: {file_path}")

        try:
            # Read CSV with custom separator
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
                # Combine values for target_id and score_rank
                data_dict[query_id]['target_id'] = data_dict[query_id]['target_id'] + "; " + str(row['target_id'])
                data_dict[query_id]['score_rank'] = str(data_dict[query_id]['score_rank']) + "; " + str(row['score_rank'])
                # Append the sample to the list
                if sample not in data_dict[query_id]['sample']:
                    data_dict[query_id]['sample'].append(sample)
            else:
                # Create a new entry in the dictionary
                data_dict[query_id] = {'target_id': row['target_id'], 'score_rank': row['score_rank'], 'sample': [sample]}
                
                # Extract additional columns dynamically
                additional_columns = annotation_data.columns.difference([query_id_col, 'target_id', 'score_rank'])
                for col in additional_columns:
                    # Assign correct values to each additional column
                    if col.lower() == 'family':
                        data_dict[query_id][col] = row['subfam-GenBank']
                    elif col.lower() == 'bitscore':
                        data_dict[query_id][col] = row['subfam-EC']
                    else:
                        data_dict[query_id][col] = row[col]

            logging.info(f"Processed query_id: {query_id} for sample: {sample}")

    # Create a DataFrame from the dictionary
    combined_data = pd.DataFrame.from_dict(data_dict, orient='index')
    combined_data.reset_index(inplace=True)
    
    # Rearrange the columns to match the desired order
    output_columns_order = [query_id_col, 'sample', 'target_id', 'score_rank', 'bitScore'] + list(all_columns.difference([query_id_col, 'target_id', 'score_rank', 'subfam-GenBank', 'subfam-EC']))
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

    combine_annotations(args.annotations, args.output)
