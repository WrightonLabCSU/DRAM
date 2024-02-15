import sys
import pandas as pd

def main(sample, db_name, descriptions_path, bit_score_threshold):
    # Load the MMseqs output
    mmseqs_path = f"mmseqs_out/{sample}_mmseqs_{db_name}.tsv"
    df_mmseqs = pd.read_csv(mmseqs_path, sep='\t', header=None)
    df_mmseqs.columns = ['query_id', 'subject_id', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

    # Rename and select columns
    df_mmseqs = df_mmseqs[['query_id', 'qstart', 'qend', 'subject_id', 'bitscore']]
    df_mmseqs.columns = ['query_id', 'start_position', 'end_position', 'merops_id', 'merops_bitScore']

    # Load the descriptions file
    if descriptions_path != "NULL":
        df_descriptions = pd.read_csv(descriptions_path, sep='\t', header=None)
        # Assuming the first column contains the gene id value
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
