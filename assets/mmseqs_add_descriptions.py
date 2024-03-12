import sys
import pandas as pd

def assign_rank(row, a_rank, b_rank, db_name_bit_score):
    """Assign canthyd rank based on bit score and provided thresholds."""
    if pd.isna(row[db_name_bit_score]):
        return None
    elif row[db_name_bit_score] >= a_rank:
        return 'A'
    elif row[db_name_bit_score] >= b_rank:
        return 'B'
    else:
        return None

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
        if first_line != "NULL":
            df_descriptions = pd.read_csv(descriptions_path, sep='\t')
            
            # Check if A_rank and B_rank columns are present
            if 'A_rank' in df_descriptions.columns and 'B_rank' in df_descriptions.columns:
                # Add database name prefix to 'bitScore' to match merged DataFrame
                db_name_bit_score = f"{db_name}_bitScore"
                
                # Ensure descriptions columns are prefixed with db_name
                df_descriptions = df_descriptions.rename(columns=lambda x: f"{db_name}_{x}" if x not in ['A_rank', 'B_rank'] else x)
                
                # Merge descriptions DataFrame
                df_merged = pd.merge(df_merged, df_descriptions, on=f'{db_name}_id', how='left')
                
                # Assign ranks based on A_rank and B_rank thresholds
                df_merged[f'{db_name}_rank'] = df_merged.apply(lambda row: assign_rank(row, row['A_rank'], row['B_rank'], db_name_bit_score), axis=1)
                
                # Clean up: Drop A_rank and B_rank columns if you don't need them anymore
                df_merged.drop(columns=['A_rank', 'B_rank'], inplace=True)
            else:
                # Handle case where descriptions do not contain A_rank and B_rank
                print("Descriptions do not contain 'A_rank' and 'B_rank', skipping rank assignment.")
    

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
