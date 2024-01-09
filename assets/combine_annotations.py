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

        # Identify columns based on suffixes
        target_id_col = identify_columns(annotation_data.columns, "_id")
        score_rank_col = identify_columns(annotation_data.columns, "_score_rank")
        bitScore_col = identify_columns(annotation_data.columns, "_bitScore")

        for index, row in annotation_data.iterrows():
            try:
                query_id = row['query_id']
            except KeyError as e:
                logging.error(f"KeyError: {e} in row: {row}")
                continue

            # Check if query_id already exists in the dictionary
            if query_id in data_dict:
                # Keep the duplicate with the highest bitScore
                if row[bitScore_col] > data_dict[query_id][bitScore_col]:
                    data_dict[query_id] = row.to_dict()
            else:
                data_dict[query_id] = row.to_dict()

        logging.info(f"Processed annotation file: {file_path}")

    # Create a DataFrame from the dictionary
    combined_data = pd.DataFrame.from_dict(data_dict, orient='index')
    combined_data.reset_index(inplace=True)
    
    # Reorder columns alphabetically
    combined_data = combined_data.reindex(sorted(combined_data.columns), axis=1)

    # Remove duplicate samples within the same row and separate them with a semicolon
    combined_data['sample'] = combined_data['sample'].apply(lambda x: "; ".join(list(set(x))))

    # Select only the required columns in the output
    output_columns = ['query_id', 'sample'] + [col for col in combined_data.columns if col != 'query_id']
    combined_data = combined_data[output_columns]

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
