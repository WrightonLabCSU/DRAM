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
    # Base columns are defined and should always appear first
    base_columns = ['query_id', 'sample', 'start_position', 'stop_position', 'strandedness', 'rank']
    # Identify all columns that start with 'kegg_'
    kegg_columns = [col for col in df.columns if col.startswith('kegg_')]
    # Identify other columns that are not in base_columns or kegg_columns
    other_columns = [col for col in df.columns if col not in base_columns and col not in kegg_columns]
    # The final order starts with base_columns, followed by kegg_columns, then the rest
    final_columns_order = base_columns + kegg_columns + other_columns
    # Reorder the DataFrame according to the final column order and return
    return df[final_columns_order]

def combine_annotations(annotation_files, output_file, threads):
    samples_and_paths = [(annotation_files[i].strip('[], '), annotation_files[i + 1].strip('[], ')) for i in range(0, len(annotation_files), 2)]
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(read_and_preprocess, sample, path) for sample, path in samples_and_paths]
        combined_data = pd.concat([future.result() for future in as_completed(futures)], ignore_index=True)

    combined_data = combined_data.drop_duplicates(subset=['query_id', 'start_position', 'stop_position'])
    combined_data = convert_bit_scores_to_numeric(combined_data)
    combined_data['rank'] = combined_data.apply(assign_rank, axis=1)
    combined_data = organize_columns(combined_data)
    combined_data.sort_values(by=['sample', 'query_id'], ascending=True, inplace=True)

    # Extract unique gene identifier from 'query_id'
    combined_data['gene_identifier'] = combined_data['query_id'].apply(lambda x: "_".join(x.split('_')[:-1]))
    # Group by 'sample' and 'gene_identifier' and enumerate entries
    combined_data['gene_number'] = combined_data.groupby(['sample', 'gene_identifier']).cumcount() + 1
    # Drop the temporary 'gene_identifier' column
    combined_data.drop(columns=['gene_identifier'], inplace=True)

    # Ensure "gene_number" comes after "rank"
    # Define the desired column order
    col_order = combined_data.columns.tolist()
    # Move 'gene_number' to the position right after 'rank'
    rank_index = col_order.index('rank')
    col_order.insert(rank_index + 1, col_order.pop(col_order.index('gene_number')))
    # Reorder the DataFrame according to the new column order
    combined_data = combined_data[col_order]

    combined_data.to_csv(output_file, index=False, sep='\t')
    logging.info(f"Combined annotations saved to {output_file}, with ranks assigned, sorted by 'sample' and 'query_id', and 'gene_number' correctly positioned.")

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
