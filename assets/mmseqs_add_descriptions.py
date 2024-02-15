import pandas as pd

def main(sample, db_name, descriptions_path, bit_score_threshold):
    # Read the output TSV file
    tsv_path = f"mmseqs_out/{sample}_mmseqs_{db_name}.tsv"
    df_mmseqs = pd.read_csv(tsv_path, sep='\t', header=None)

    # Rename columns
    df_mmseqs.columns = ['query_id', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

    # Select and rename desired columns
    df_mmseqs = df_mmseqs[['query_id', 'qstart', 'qend', 'sseqid', 'bitscore']]
    df_mmseqs.columns = ['query_id', 'start_position', 'end_position', 'merops_id', 'merops_bitScore']

    # Merge with database descriptions if available
    if descriptions_path != "NULL":
        df_descriptions = pd.read_csv(descriptions_path, sep='\t')
        # Rename columns with db_name prefix
        df_descriptions = df_descriptions.add_prefix(f"{db_name}_")
        # Merge on query_id
        df_merged = pd.merge(df_mmseqs, df_descriptions, on='query_id', how='left')
        # Output merged dataframe to CSV
        output_path = f"mmseqs_out/{sample}_mmseqs_{db_name}_formatted.csv"
        df_merged.to_csv(output_path, index=False)
    else:
        # Output parsed dataframe to CSV
        output_path = f"mmseqs_out/{sample}_mmseqs_{db_name}_formatted.csv"
        df_mmseqs.to_csv(output_path, index=False)

if __name__ == "__main__":
    import sys
    sample = sys.argv[1]
    db_name = sys.argv[2]
    descriptions_path = sys.argv[3]
    bit_score_threshold = float(sys.argv[4])
    main(sample, db_name, descriptions_path, bit_score_threshold)
