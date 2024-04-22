import pandas as pd
import sys

def extract_query_ids(tsv_path, ko_list):
    df = pd.read_csv(tsv_path, sep='\t')
    
    # Split the ko_list string into terms, separating on semicolons
    ko_terms = ko_list.split(';')
    
    # Identify columns for gene IDs, EC numbers, and descriptions
    id_columns = [col for col in df.columns if col.endswith('_id') and col != 'query_id']
    ec_columns = [col for col in df.columns if col.endswith('_EC')]
    desc_columns = [col for col in df.columns if col.endswith('_description')]
    
    # Create filters for ID, EC, and description columns
    id_mask = df[id_columns].apply(lambda x: x.isin(ko_terms)).any(axis=1)
    ec_mask = df[ec_columns].apply(lambda x: x.isin(ko_terms)).any(axis=1)
    desc_mask = df[desc_columns].apply(lambda x: x.astype(str).apply(lambda y: any(term in y for term in ko_terms))).any(axis=1)
    
    # Combine all masks to filter the DataFrame
    combined_mask = id_mask | ec_mask | desc_mask
    filtered_df = df[combined_mask]
    
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
