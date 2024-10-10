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

def organize_columns(df, special_columns=None):
    if special_columns is None:
        special_columns = []
    base_columns = ['query_id', 'sample', 'start_position', 'stop_position', 'strandedness', 'rank', 'gene_number']
    base_columns = [col for col in base_columns if col in df.columns]
    
    kegg_columns = sorted([col for col in df.columns if col.startswith('kegg_')], key=lambda x: (x != 'kegg_id', x))
    other_columns = [col for col in df.columns if col not in base_columns + kegg_columns + special_columns]
    
    db_prefixes = set(col.split('_')[0] for col in other_columns)
    sorted_other_columns = []
    for prefix in db_prefixes:
        prefixed_columns = sorted([col for col in other_columns if col.startswith(prefix + '_')], key=lambda x: (x != f"{prefix}_id", x))
        sorted_other_columns.extend(prefixed_columns)
    
    final_columns_order = base_columns + kegg_columns + sorted_other_columns + special_columns
    return df[final_columns_order]


def combine_annotations(annotation_files, output_file, threads):
    samples_and_paths = [(annotation_files[i].strip('[], '), annotation_files[i + 1].strip('[], ')) for i in range(0, len(annotation_files), 2)]
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(read_and_preprocess, sample, path) for sample, path in samples_and_paths]
        data_frames = [future.result() for future in as_completed(futures)]
    
    combined_data = pd.concat(data_frames, ignore_index=True)
    combined_data = convert_bit_scores_to_numeric(combined_data)
    
    aggregation_functions = {col: 'first' for col in combined_data.columns if col not in ['query_id', 'sample']}
    for col in ['Completeness', 'Contamination', 'taxonomy']:
        if col in combined_data.columns:
            aggregation_functions[col] = 'max'
    
    combined_data = combined_data.groupby(['query_id', 'sample'], as_index=False).agg(aggregation_functions)
    # After aggregating data
    combined_data['rank'] = combined_data.apply(assign_rank, axis=1)

    # Correctly extract the base part of 'query_id'
    combined_data['base_query_id'] = combined_data['query_id'].str.rsplit('_', n=1).str[0]
    # Recalculate 'gene_number' with corrected grouping
    combined_data['gene_number'] = combined_data.groupby(['sample', 'base_query_id']).cumcount() + 1

    # Continue with organizing columns and saving the DataFrame
    special_columns = ['Completeness', 'Contamination', 'taxonomy']
    columns_to_exclude = [col for col in special_columns if col in combined_data.columns]
    combined_data = organize_columns(combined_data, special_columns=columns_to_exclude)

    combined_data.drop(columns=['base_query_id'], inplace=True)  # Cleanup

    combined_data.to_csv(output_file, index=False, sep='\t')
    logging.info(f"Combined annotations saved to {output_file}, with corrected gene numbers.")



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
