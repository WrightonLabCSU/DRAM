import argparse
import pandas as pd
import logging

# Configure the logger
logging.basicConfig(filename="logs/combine_annotations.log", level=logging.INFO, format='%(levelname)s: %(message)s')

def combine_annotations(annotation_files, output_file):
    # Create an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()

    # Create a dictionary to store values for each unique query_id
    data_dict = {}

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

        for index, row in annotation_data.iterrows():
            try:
                query_id = row['query_id']
            except KeyError as e:
                logging.error(f"KeyError: {e} in row: {row}")
                continue

            # Check if query_id already exists in the dictionary
            if query_id in data_dict:
                # Combine values for target_id and score_rank
                data_dict[query_id]['target_id'] = row.get('_id', '')
                data_dict[query_id]['score_rank'] = row.get('_score_rank', '')
                # Append the sample to the list
                if sample not in data_dict[query_id]['sample']:
                    data_dict[query_id]['sample'].append(sample)
            else:
                # Create a new entry in the dictionary
                data_dict[query_id] = {'target_id': row.get('_id', ''), 'score_rank': row.get('_score_rank', ''), 'sample': [sample]}
                
                # Retain additional columns in the output
                additional_columns = annotation_data.columns.difference(['query_id', '_id', '_score_rank'])
                for col in additional_columns:
                    data_dict[query_id][col] = row[col]

            logging.info(f"Processed query_id: {query_id} for sample: {sample}")

    # Create a DataFrame from the dictionary
    combined_data = pd.DataFrame.from_dict(data_dict, orient='index')
    combined_data.reset_index(inplace=True)
    
    # Reorder columns
    combined_data = combined_data[['query_id', 'target_id', 'sample', 'score_rank'] + sorted(list(additional_columns))]

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
