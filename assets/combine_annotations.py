import argparse
import pandas as pd
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure the logger
logging.basicConfig(filename="logs/combine_annotations.log", level=logging.INFO, format='%(levelname)s: %(message)s')

def read_and_preprocess(sample, path):
    try:
        df = pd.read_csv(path)
        df['sample'] = sample  # Add sample column
        return df
    except Exception as e:
        logging.error(f"Error loading DataFrame for sample {sample}: {str(e)}")
        return pd.DataFrame()  # Return an empty DataFrame in case of error

def assign_rank(row):
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

def organize_columns(df):
    # Base columns to always appear first
    base_columns = ['query_id', 'sample', 'start_position', 'stop_position', 'strandedness', 'rank', 'gene_number']
    
    # Extract KEGG columns and ensure they are listed first if present
    kegg_columns = [col for col in df.columns if col.startswith('kegg_')]
    kegg_columns_sorted = sorted(kegg_columns, key=lambda x: (x != 'kegg_id', x))
    
    # Identify and organize other database columns
    db_columns = [col for col in df.columns if col not in base_columns and col not in kegg_columns]
    
    # Identify unique database prefixes
    db_prefixes = set(col.split('_')[0] for col in db_columns)
    
    # Sort and organize columns for each database, ensuring <db_name>_id is the first column
    other_db_columns_sorted = []
    for prefix in sorted(db_prefixes):
        prefixed_columns = [col for col in db_columns if col.startswith(prefix + '_')]
        prefixed_columns_sorted = sorted(prefixed_columns, key=lambda x: (x != f"{prefix}_id", x))
        other_db_columns_sorted.extend(prefixed_columns_sorted)
    
    # Final order of columns
    final_columns_order = base_columns + kegg_columns_sorted + other_db_columns_sorted
    return df[final_columns_order]


def combine_annotations(annotation_files, output_file, threads):
    samples_and_paths = [(annotation_files[i].strip('[], '), annotation_files[i + 1].strip('[], ')) for i in range(0, len(annotation_files), 2)]
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(read_and_preprocess, sample, path) for sample, path in samples_and_paths]
        data_frames = [future.result() for future in as_completed(futures)]
    
    # Initialize an empty DataFrame for combined data
    combined_data = pd.DataFrame()
    # Track if special columns have been added already
    special_columns_added = {'Completeness': False, 'Contamination': False, 'taxonomy': False}
    
    for df in data_frames:
        # Check and drop special columns if they've been added already
        for col in ['Completeness', 'Contamination', 'taxonomy']:
            if col in df.columns and special_columns_added[col]:
                df.drop(columns=[col], inplace=True)
            elif col in df.columns:
                special_columns_added[col] = True
                
        combined_data = pd.concat([combined_data, df], ignore_index=True)
    
    combined_data = combined_data.drop_duplicates(subset=['query_id', 'start_position', 'stop_position', 'sample'])
    combined_data = convert_bit_scores_to_numeric(combined_data)
    combined_data['rank'] = combined_data.apply(assign_rank, axis=1)
    
    combined_data = organize_columns(combined_data)
    combined_data['gene_number'] = combined_data.groupby(['sample', 'query_id']).cumcount() + 1
    
    combined_data.to_csv(output_file, index=False, sep='\t')
    logging.info(f"Combined annotations saved to {output_file}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine annotation files with ranks and avoid duplicating specific columns.")
    parser.add_argument("--annotations", nargs='+', help="List of annotation files and sample names, alternating.")
    parser.add_argument("--threads", help="Number of threads for parallel processing", type=int, default=4)
    parser.add_argument("--output", help="Output file path for the combined annotations.")
    args = parser.parse_args()

    if args.annotations and args.output:
        combine_annotations(args.annotations, args.output, args.threads)
    else:
        logging.error("Missing required arguments. Use --help for usage information.")
