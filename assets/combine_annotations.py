import argparse
import pandas as pd
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure the logger
logging.basicConfig(filename="logs/combine_annotations.log", level=logging.INFO, format='%(levelname)s: %(message)s')

def read_and_preprocess(sample_path):
    sample, path = sample_path
    try:
        df = pd.read_csv(path)
        df['sample'] = sample  # Add sample column
        # Pre-process data here if needed
        return df
    except Exception as e:
        logging.error(f"Error loading DataFrame for sample {sample}: {str(e)}")
        return pd.DataFrame()  # Return an empty DataFrame in case of error

def assign_rank(row):
    # Rank assignment logic remains unchanged
    rank = 'E'
    if row.get('kegg_bitScore', 0) > 350:
        rank = 'A'
    elif row.get('uniref_bitScore', 0) > 350:
        rank = 'B'
    elif row.get('kegg_bitScore', 0) > 60 or row.get('uniref_bitScore', 0) > 60:
        rank = 'C'
    elif any(row.get(f"{db}_bitScore", 0) > 60 for db in ['pfam', 'dbcan', 'merops']):
        rank = 'D'
    return rank

def convert_bit_scores_to_numeric(df):
    for col in df.columns:
        if "_bitScore" in col:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def combine_annotations(annotation_files, output_file):
    samples_and_paths = [(annotation_files[i], annotation_files[i + 1]) for i in range(0, len(annotation_files), 2)]
    
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(read_and_preprocess, sp) for sp in samples_and_paths]
        combined_data = pd.concat([f.result() for f in as_completed(futures)], ignore_index=True)

    # Post-processing
    combined_data.drop_duplicates(subset=['query_id', 'start_position', 'stop_position'], inplace=True)
    combined_data = convert_bit_scores_to_numeric(combined_data)
    combined_data['rank'] = combined_data.apply(assign_rank, axis=1)

    # Organize columns
    base_columns = ['query_id', 'sample', 'start_position', 'stop_position', 'strandedness', 'rank']
    other_columns = [col for col in combined_data.columns if col not in base_columns]
    combined_data = combined_data[base_columns + other_columns]
    
    # Sort and save
    combined_data.sort_values(by='query_id', ascending=True, inplace=True)
    combined_data.to_csv(output_file, index=False, sep='\t')
    logging.info(f"Combined annotations saved to {output_file}, with ranks assigned and sorted by 'query_id'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine annotation files with ranks and avoid duplicating specific columns.")
    parser.add_argument("--annotations", nargs='+', help="List of annotation files and sample names, alternating.")
    parser.add_argument("--threads", help="Number of threads for parallel processing", type=int, default=4)
    parser.add_argument("--output", help="Output file path for the combined annotations.")
    args = parser.parse_args()

    if args.annotations and args.output:
        combine_annotations(args.annotations, args.output)
    else:
        logging.error("Missing required arguments. Use --help for usage information.")
