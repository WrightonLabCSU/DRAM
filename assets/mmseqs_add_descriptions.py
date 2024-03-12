import sys
import pandas as pd

def assign_rank(row, a_rank, b_rank, db_name_bit_score):
    """Assign rank based on bit score and provided thresholds."""
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

    mmseqs_path = f"mmseqs_out/{sample}_mmseqs_{db_name}.tsv"
    print(f"Loading MMseqs output from {mmseqs_path}")
    
    df_gene_locs = pd.read_csv(gene_locs_path, sep='\t', header=None, names=['query_id', 'start_position', 'stop_position'])
    print("Gene locations loaded. Sample rows:")
    print(df_gene_locs.head())
    
    df_mmseqs = pd.read_csv(mmseqs_path, sep='\t', header=None, usecols=[0, 1, 11])
    df_mmseqs.columns = ['query_id', f'{db_name}_id', f'{db_name}_bitScore']
    print("MMseqs data loaded and columns renamed. Sample rows:")
    print(df_mmseqs.head())
    
    df_merged = pd.merge(df_mmseqs, df_gene_locs, on='query_id', how='left')
    print("Merged MMseqs output with gene locations. Sample rows:")
    print(df_merged.head())

    if descriptions_path != "NULL":
        df_descriptions = pd.read_csv(descriptions_path, sep='\t')

        # Rename 'definition' to 'description' and prefix other columns, excluding 'A_rank' and 'B_rank'
        df_descriptions.rename(columns={'definition': 'description'}, inplace=True)
        df_descriptions = df_descriptions.rename(columns={col: f"{db_name}_{col}" for col in df_descriptions.columns if col not in ['A_rank', 'B_rank', f'{db_name}_id']})

        # Ensure the first column is renamed to match db_name_id for merging
        first_column = df_descriptions.columns[0]
        df_descriptions.rename(columns={first_column: f'{db_name}_id'}, inplace=True)

        df_merged = pd.merge(df_merged, df_descriptions, on=f'{db_name}_id', how='left')

        db_name_bit_score = f"{db_name}_bitScore"
        df_merged[f'{db_name}_rank'] = df_merged.apply(lambda row: assign_rank(row, row[f'{db_name}_A_rank'], row[f'{db_name}_B_rank'], db_name_bit_score), axis=1)
        
        # Drop A_rank and B_rank columns if present
        df_merged.drop(columns=[f'{db_name}_A_rank', f'{db_name}_B_rank'], inplace=True, errors='ignore')

    output_path = f"mmseqs_out/{sample}_mmseqs_{db_name}_formatted.csv"
    df_merged.to_csv(output_path, index=False)
    print("Merged DataFrame saved to", output_path)

if __name__ == "__main__":
    sample = sys.argv[1]
    db_name = sys.argv[2]
    descriptions_path = sys.argv[3]
    bit_score_threshold = float(sys.argv[4])
    gene_locs_path = sys.argv[5]

    main(sample, db_name, descriptions_path, bit_score_threshold, gene_locs_path)
