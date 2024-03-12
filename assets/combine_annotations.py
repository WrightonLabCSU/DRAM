import argparse
import pandas as pd
import logging

# Configure the logger
logging.basicConfig(filename="logs/combine_annotations.log", level=logging.INFO, format='%(levelname)s: %(message)s')

def extract_samples_and_paths(annotation_files):
    samples_and_paths = []
    for i in range(0, len(annotation_files), 2):
        sample = annotation_files[i].strip('[], ')
        path = annotation_files[i + 1].strip('[], ')
        samples_and_paths.append((sample, path))
    return samples_and_paths

def organize_columns(df):
    # Define base columns and ensure they are at the beginning
    base_columns = ['query_id', 'sample', 'start_position', 'stop_position', 'strandedness']
    base_columns_present = [col for col in base_columns if col in df.columns]

    # Find and prioritize KEGG columns if present
    kegg_columns = [col for col in df.columns if col.startswith('kegg_')]
    
    # Find other database columns
    other_db_columns = [col for col in df.columns if col not in base_columns_present + kegg_columns]
    
    # The final column order starts with base columns, followed by KEGG columns, then the rest
    final_columns_order = base_columns_present + kegg_columns + other_db_columns
    
    return df[final_columns_order]

def combine_annotations(annotation_files, output_file):
    # Extract samples and paths from the input annotation_files
    samples_and_paths = extract_samples_and_paths(annotation_files)

    # Create an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()

    # Process the input annotation files
    for sample, path in samples_and_paths:
        # Load each annotation file into a DataFrame
        try:
            annotation_df = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            logging.warning(f"Empty DataFrame for sample: {sample}, skipping.")
            continue
        except Exception as e:
            logging.error(f"Error loading DataFrame for sample {sample}: {str(e)}")
            continue

        # Add the 'sample' column
        annotation_df.insert(0, 'sample', sample)

        # Combine data
        combined_data = pd.concat([combined_data, annotation_df], ignore_index=True, sort=False)

    # Modify grouping to preserve unique combinations of 'query_id', 'start_position', and 'stop_position'
    combined_data = combined_data.drop_duplicates(subset=['query_id', 'start_position', 'stop_position'])

    # Reorder the columns according to the new requirements
    combined_data = organize_columns(combined_data)

    # Save the combined DataFrame to the output file
    combined_data.to_csv(output_file, index=False, sep='\t')

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Combine annotation files preserving unique combinations of 'query_id', 'start_position', 'stop_position', and including 'strandedness'. Prioritize 'kegg' database columns if present.")
    parser.add_argument("--annotations", nargs='+', help="List of annotation files and sample names.")
    parser.add_argument("--output", help="Output file path for the combined annotations.")
    args = parser.parse_args()

    # Combine annotations
    if args.annotations and args.output:
        combine_annotations(args.annotations, args.output)
        logging.info(f"Combined annotations saved to {args.output}")
    else:
        logging.error("Missing required arguments. Use --help for usage information.")
