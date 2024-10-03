import pandas as pd
import sys

def extract_query_ids(tsv_path, ko_terms):
    df = pd.read_csv(tsv_path, sep='\t')
    
    # Split KO terms on newline to handle multiple entries correctly
    ko_terms = ko_terms.strip().split('\n')
    print(f"Searching for KOs: {ko_terms}")

    # Define columns to search for KO terms, excluding 'query_id'
    search_columns = [col for col in df.columns if (col.endswith('_id') or col.endswith('_description') or col.endswith('_EC')) and col != 'query_id']
    
    # Debugging: Output the columns being searched
    print(f"Searching in columns: {search_columns}")

    # Create a filter mask for any row containing any of the KO terms in the specified columns
    mask = df[search_columns].apply(lambda x: x.astype(str).str.contains('|'.join(ko_terms), case=False, na=False)).any(axis=1)
    filtered_df = df[mask]

    print(f"Found {len(filtered_df)} matching entries.")

    return filtered_df[['sample', 'query_id']]

def main():
    if len(sys.argv) != 4:
        print("Usage: python parse_annotations.py <tsv_path> <ko_list> <output_file>")
        sys.exit(1)

    tsv_path, ko_list_path, output_file = sys.argv[1:]
    with open(ko_list_path, 'r') as f:
        ko_list = f.read()

    results_df = extract_query_ids(tsv_path, ko_list)

    # Write results to output file
    results_df.to_csv(output_file, sep='\t', index=False, header=False)

    print(f"Extracted {len(results_df)} entries written to {output_file}")

if __name__ == '__main__':
    main()
