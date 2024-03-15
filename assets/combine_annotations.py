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
    base_columns = ['query_id', 'sample', 'start_position', 'stop_position', 'strandedness', 'rank']
    kegg_columns = [col for col in df.columns if col.startswith('kegg_')]
    other_columns = [col for col in df.columns if col not in base_columns and col not in kegg_columns]
    final_columns_order = base_columns + kegg_columns + other_columns
    return df[final_columns_order]

def combine_annotations(annotation_files, output_file, threads):
    samples_and_paths = [(annotation_files[i].strip('[], '), annotation_files[i + 1].strip('[], ')) for i in range(0, len(annotation_files), 2)]
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(read_and_preprocess, sample, path) for sample, path in samples_and_paths]
        combined_data = pd.concat([future.result() for future in as_completed(futures)], ignore_index=True)

    combined_data = combined_data.drop_duplicates(subset=['query_id', 'start_position', 'stop_position'])
    combined_data = convert_bit_scores_to_numeric(combined_data)
    combined_data['rank'] = combined_data.apply(assign_rank, axis=1)
    
    # Create a new identifier for sorting that combines the non-numeric part of query_id and its numeric suffix
    combined_data['sorting_id'] = combined_data['query_id'].str.extract(r'(.+)_([0-9]+)$', expand=True).apply(lambda x: f"{x[0]}_{int(x[1]):05d}", axis=1)
    
    combined_data = organize_columns(combined_data)
    # Sort data by 'sample' and 'sorting_id' to correctly assign 'gene_number'
    combined_data.sort_values(by=['sample', 'sorting_id'], ascending=True, inplace=True)

    # Group by 'sample' and original 'query_id' base (without the numeric part) and enumerate entries based on their sorted order
    combined_data['base_query_id'] = combined_data['query_id'].apply(lambda x: "_".join(x.split('_')[:-1]))
    combined_data['gene_number'] = combined_data.groupby(['sample', 'base_query_id']).cumcount() + 1
    
    # Ensure "gene_number" comes after "rank"
    col_order = combined_data.columns.tolist()
    rank_index = col_order.index('rank')
    col_order.insert(rank_index + 1, col_order.pop(col_order.index('gene_number')))
    combined_data = combined_data[col_order]

    combined_data.drop(columns=['sorting_id', 'base_query_id'], inplace=True)  # Drop the temporary 'sorting_id' and 'base_query_id' columns

    combined_data.to_csv(output_file, index=False, sep='\t')
    logging.info(f"Combined annotations saved to {output_file}, with ranks assigned, sorted by 'sample' and 'query_id' numerically, and 'gene_number' correctly positioned.")

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
