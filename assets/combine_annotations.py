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
    # Dynamically check for database hits and assign ranks based on bitScore criteria
    databases = ["kegg", "uniref", "pfam", "dbcan", "merops", "vogdb"]
    db_scores = {db: row.get(f"{db}_bitScore", None) for db in databases}
    
    # Rank logic
    if db_scores["kegg"] is not None and db_scores["kegg"] > 350:
        return 'A'
    elif db_scores["uniref"] is not None and db_scores["uniref"] > 350:
        return 'B'
    elif (db_scores["kegg"] is not None and db_scores["kegg"] > 60) or (db_scores["uniref"] is not None and db_scores["uniref"] > 60):
        return 'C'
    elif any(db in row for db in ["pfam_id", "dbcan_id", "merops_id"]) and all(score <= 60 for score in db_scores.values() if score is not None):
        return 'D'
    elif db_scores["vogdb"] is not None or all(score is None or score < 60 for score in db_scores.values()):
        return 'E'
    return 'E'  # Default to rank E if no other conditions are met

def organize_columns(df):
    base_columns = ['query_id', 'sample', 'start_position', 'stop_position', 'strandedness', 'rank']
    kegg_columns = [col for col in df.columns if col.startswith('kegg_')]
    other_db_columns = [col for col in df.columns if col not in base_columns + kegg_columns]
    
    final_columns_order = base_columns + kegg_columns + other_db_columns
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

    combined_data = combined_data.drop_duplicates(subset=['query_id', 'start_position', 'stop_position'])
    combined_data['rank'] = combined_data.apply(assign_rank, axis=1)
    combined_data = organize_columns(combined_data)
    combined_data = combined_data.sort_values(by='query_id', ascending=True)

    combined_data.to_csv(output_file, index=False, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine annotation files with ranking based on annotation confidence and sort by 'query_id'.")
    parser.add_argument("--annotations", nargs='+', help="List of annotation files and sample names.")
    parser.add_argument("--output", help="Output file path for the combined annotations.")
    args = parser.parse_args()

    if args.annotations and args.output:
        combine_annotations(args.annotations, args.output)
        logging.info(f"Combined annotations saved to {args.output}, with ranks assigned and sorted by 'query_id'.")
    else:
        logging.error("Missing required arguments. Use --help for usage information.")
