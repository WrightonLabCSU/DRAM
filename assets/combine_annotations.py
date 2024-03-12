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

def assign_rank(row):
    # Initial rank is 'E'
    rank = 'E'
    
    # Logic to assign ranks based on bit scores from various databases
    if row.get('kegg_bitScore', 0) > 350:
        rank = 'A'
    elif row.get('uniref_bitScore', 0) > 350:
        rank = 'B'
    elif row.get('kegg_bitScore', 0) > 60 or row.get('uniref_bitScore', 0) > 60:
        rank = 'C'
    elif any(row.get(f"{db}_bitScore", 0) > 60 for db in ['pfam', 'dbcan', 'merops']):
        rank = 'D'
    # No need for an explicit check for rank 'E', as it's the default

    # Optional: Log the decision for auditing or debugging
    # logging.debug(f"Assigned rank {rank} for query_id {row.get('query_id', 'N/A')} with kegg_bitScore {row.get('kegg_bitScore', 'N/A')}")

    return rank

def convert_bit_scores_to_numeric(df):
    for col in df.columns:
        if "_bitScore" in col:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def organize_columns(df):
    # Ensure the 'rank' column is included after 'strandedness'
    base_columns = ['query_id', 'sample', 'start_position', 'stop_position', 'strandedness', 'rank']
    other_columns = [col for col in df.columns if col not in base_columns]
    final_columns_order = base_columns + other_columns
    return df[final_columns_order]

def combine_annotations(annotation_files, output_file):
    samples_and_paths = extract_samples_and_paths(annotation_files)
    combined_data = pd.DataFrame()

    for sample, path in samples_and_paths:
        try:
            annotation_df = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            logging.warning(f"Empty DataFrame for sample: {sample}, skipping.")
            continue
        except Exception as e:
            logging.error(f"Error loading DataFrame for sample {sample}: {str(e)}")
            continue

        annotation_df.insert(0, 'sample', sample)
        combined_data = pd.concat([combined_data, annotation_df], ignore_index=True, sort=False)

    # Drop duplicates
    combined_data = combined_data.drop_duplicates(subset=['query_id', 'start_position', 'stop_position'])

    # Convert bit scores to numeric for accurate comparisons
    combined_data = convert_bit_scores_to_numeric(combined_data)

    # Assign ranks based on defined criteria after ensuring numeric conversion
    combined_data['rank'] = combined_data.apply(assign_rank, axis=1)

    # Organize and sort columns, including the newly added 'rank' column
    combined_data = organize_columns(combined_data)
    combined_data = combined_data.sort_values(by='query_id', ascending=True)

    combined_data.to_csv(output_file, index=False, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine annotation files with ranks and avoid duplicating specific columns.")
    parser.add_argument("--annotations", nargs='+', help="List of annotation files and sample names.")
    parser.add_argument("--output", help="Output file path for the combined annotations.")
    args = parser.parse_args()

    if args.annotations and args.output:
        combine_annotations(args.annotations, args.output)
        logging.info(f"Combined annotations saved to {args.output}, with ranks assigned and sorted by 'query_id'.")
    else:
        logging.error("Missing required arguments. Use --help for usage information.")
