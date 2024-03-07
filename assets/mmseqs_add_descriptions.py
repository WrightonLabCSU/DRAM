import sys
import pandas as pd

def main(sample, db_name, descriptions_path, bit_score_threshold, gene_locs_path):
    # Load the MMseqs output
    mmseqs_path = f"mmseqs_out/{sample}_mmseqs_{db_name}.tsv"
    
    # Load gene locations
    df_gene_locs = pd.read_csv(gene_locs_path, sep='\t', header=None, names=['query_id', 'start_position', 'stop_position'])
    
    # Read the MMseqs file with the necessary columns
    df_mmseqs = pd.read_csv(mmseqs_path, sep='\t', header=None, usecols=[0, 1, 11])
    
    # Rename the columns accordingly
    df_mmseqs.columns = ['query_id', f'{db_name}_id', f'{db_name}_bitScore']
    
    # Merge MMseqs output with gene locations based on 'query_id'
    df_merged = pd.merge(df_mmseqs, df_gene_locs, on='query_id', how='left')
    
    # Load the descriptions file if it's provided
    if descriptions_path != "NULL":
        df_descriptions = pd.read_csv(descriptions_path, sep='\t')
        
        # Add database name prefix to added columns
        df_descriptions.columns = [f"{db_name}_{col}" if col != 'blast_name' else col for col in df_descriptions.columns]
        
        print("Columns in df_merged:", df_merged.columns)
        print("Columns in df_descriptions:", df_descriptions.columns)
        
        # Merge the DataFrames on the database ID
        df_merged = pd.merge(df_merged, df_descriptions, left_on=f"{db_name}_id", right_on=f"{db_name}_blast_name", how='left')
        
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
