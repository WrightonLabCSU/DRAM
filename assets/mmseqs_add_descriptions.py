import sys
import pandas as pd

def main(sample, db_name, descriptions_path, bit_score_threshold, gene_locs_path):
    print(f"Starting processing for sample: {sample}, database: {db_name}")

    # Load the MMseqs output
    mmseqs_path = f"mmseqs_out/{sample}_mmseqs_{db_name}.tsv"
    print(f"Loading MMseqs output from {mmseqs_path}")
    
    # Load gene locations
    df_gene_locs = pd.read_csv(gene_locs_path, sep='\t', header=None, names=['query_id', 'start_position', 'stop_position'])
    print("Gene locations loaded. Sample rows:")
    print(df_gene_locs.head())
    
    # Read the MMseqs file with the necessary columns
    df_mmseqs = pd.read_csv(mmseqs_path, sep='\t', header=None, usecols=[0, 1, 11])
    df_mmseqs.columns = ['query_id', f'{db_name}_id', f'{db_name}_bitScore']
    print("MMseqs data loaded and columns renamed. Sample rows:")
    print(df_mmseqs.head())
    
    # Merge MMseqs output with gene locations based on 'query_id'
    df_merged = pd.merge(df_mmseqs, df_gene_locs, on='query_id', how='left')
    print("Merged MMseqs output with gene locations. Sample rows:")
    print(df_merged.head())

    # Load the descriptions file if it's provided and check its contents
    if descriptions_path != "NULL":
        with open(descriptions_path, 'r') as file:
            first_line = file.readline().strip()
        # Check if the first line of the file is "NULL"
        if first_line != "NULL":
            df_descriptions = pd.read_csv(descriptions_path, sep='\t')
            
            # Add database name prefix to added columns
            df_descriptions.columns = [f"{db_name}_{col}" for col in df_descriptions.columns]
            
            # Print statements for debugging
            print("Columns in df_merged:", df_merged.columns)
            print("Columns in df_descriptions:", df_descriptions.columns)
            
            # Merge the DataFrames on the query_id and the index of descriptions
            df_merged = pd.merge(df_merged, df_descriptions, left_index=True, right_index=True, how='left')
        else:
            print("Descriptions file content is 'NULL', skipping loading of descriptions.")
    

    # Save the merged DataFrame to CSV
    output_path = f"mmseqs_out/{sample}_mmseqs_{db_name}_formatted.csv"
    df_merged.to_csv(output_path, index=False)
    print("Merged DataFrame saved to", output_path)

if __name__ == "__main__":
    sample = sys.argv[1]
    db_name = sys.argv[2]
    descriptions_path = sys.argv[3]
    bit_score_threshold = float(sys.argv[4])  # Ensure threshold is a float
    gene_locs_path = sys.argv[5]  # Path to the gene locations file

    main(sample, db_name, descriptions_path, bit_score_threshold, gene_locs_path)
