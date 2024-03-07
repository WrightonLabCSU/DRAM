import argparse
import pandas as pd
import logging
import os

# Configure the logger
logging.basicConfig(filename="logs/combine_annotations.log", level=logging.INFO, format='%(levelname)s: %(message)s')

def extract_samples_and_paths(annotation_files):
    samples_and_paths = []
    for i in range(0, len(annotation_files), 2):
        sample = annotation_files[i].strip('[], ')
        path = annotation_files[i + 1].strip('[], ')
        samples_and_paths.append((sample, path))
    return samples_and_paths

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
        annotation_df.insert(1, 'sample', sample)

        # Concatenate the annotation DataFrame with the combined data
        combined_data = pd.concat([combined_data, annotation_df], ignore_index=True)

    # Combine annotations based on 'query_id', 'start_position', and 'stop_position'
    combined_data = combined_data.groupby(['query_id', 'start_position', 'stop_position'], as_index=False).first()

    # Save the combined DataFrame to the output file
    combined_data.to_csv(output_file, index=False, sep='\t')



if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Combine annotation files preserving unique combinations of 'query_id', 'start_position', and 'stop_position'.")
    parser.add_argument("--annotations", nargs='+', help="List of annotation files and sample names.")
    parser.add_argument("--output", help="Output file path for the combined annotations.")
    args = parser.parse_args()

    # Combine annotations
    if args.annotations and args.output:
        combine_annotations(args.annotations, args.output)
        logging.info(f"Combined annotations saved to {args.output}")
    else:
        logging.error("Missing required arguments. Use --help for usage information.")
