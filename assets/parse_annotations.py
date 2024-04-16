import pandas as pd
import sys

def extract_query_ids(tsv_path, ko_list):
    df = pd.read_csv(tsv_path, sep='\t')
    
    # Create a list of gene IDs from the provided ko_list string
    ko_terms = ko_list.split(';')
    
    # Identify columns that end with "_id" but not 'query_id'
    id_columns = [col for col in df.columns if col.endswith('_id') and col != 'query_id']
    
    # Filter DataFrame based on whether any of the id_columns contain any of the ko_terms
    mask = df[id_columns].apply(lambda x: x.isin(ko_terms)).any(axis=1)
    filtered_df = df[mask]
    
    return filtered_df[['sample', 'query_id']]

def main():
    if len(sys.argv) != 4:
        print("Usage: python parse_annotations.py <tsv_path> <ko_list> <output_file>")
        sys.exit(1)

    tsv_path, ko_list, output_file = sys.argv[1:]
    results_df = extract_query_ids(tsv_path, ko_list)

    # Write results to output file
    results_df.to_csv(output_file, sep='\t', index=False, header=False)

    print(f"Extracted {len(results_df)} entries written to {output_file}")

if __name__ == '__main__':
    main()
