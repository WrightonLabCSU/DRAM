import sys
import pandas as pd

def main(sample, db_name, descriptions_path, bit_score_threshold):
    # Load the MMseqs output
    mmseqs_path = f"mmseqs_out/{sample}_mmseqs_{db_name}.tsv"
    
    # Read the first row to determine column positions
    with open(mmseqs_path, 'r') as f:
        first_row = f.readline().strip().split('\t')

    # Determine the column positions
    query_id_index = 0
    gene_id_index = 1
    start_position_index = 6
    end_position_index = 7
    bitscore_index = 11

    # Read the file with correct column positions
    df_mmseqs = pd.read_csv(mmseqs_path, sep='\t', header=None, usecols=[query_id_index, gene_id_index, start_position_index, end_position_index, bitscore_index])
    
    # Rename the columns accordingly
    df_mmseqs.columns = ['query_id', f'{db_name}_id', 'start_position', 'end_position', f'{db_name}_bitScore']

    # Load the descriptions file
    if descriptions_path != "NULL":
        df_descriptions = pd.read_csv(descriptions_path, sep='\t', header=None)
        df_descriptions.columns = ['gene_id'] + [f'{db_name}_{i}' for i in range(1, len(df_descriptions.columns))]

        # Merge the DataFrames
        df_merged = pd.merge(df_mmseqs, df_descriptions, left_on=f"{db_name}_id", right_on='gene_id', how='left')
        
        # Save the merged DataFrame to CSV
        output_path = f"mmseqs_out/{sample}_mmseqs_{db_name}_formatted.csv"
        df_merged.to_csv(output_path, index=False)

if __name__ == "__main__":
    # Extract command-line arguments
    sample = sys.argv[1]
    db_name = sys.argv[2]
    descriptions_path = sys.argv[3]
    bit_score_threshold = float(sys.argv[4])  # Convert threshold to float

    main(sample, db_name, descriptions_path, bit_score_threshold)
