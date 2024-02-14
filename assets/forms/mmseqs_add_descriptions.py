import sys
import pandas as pd

def main(sample, db_name, descriptions_path, bit_score_threshold, cpus):
    # Define file paths
    input_path = f"mmseqs_out/{sample}_mmseqs_{db_name}.tsv"
    output_path = f"mmseqs_out/{sample}_mmseqs_{db_name}_formatted.csv"

    # Load MMseqs output into DataFrame
    mmseqs_columns = ['query_id', 'target_id', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bitscore']
    df_mmseqs = pd.read_csv(input_path, sep='\t', names=mmseqs_columns)

    # Add additional columns for start and end positions, and format target_id with db_name
    df_mmseqs['start_position'] = df_mmseqs['qstart']
    df_mmseqs['end_position'] = df_mmseqs['qend']
    df_mmseqs[f"{db_name}_id"] = df_mmseqs['target_id']
    df_mmseqs[f"{db_name}_bitScore"] = df_mmseqs['bitscore']
    df_mmseqs_final = df_mmseqs[['query_id', 'start_position', 'end_position', f"{db_name}_id", f"{db_name}_bitScore"]]

    # Check if db_descriptions is not 'NULL'
    if not pd.read_csv(descriptions_path, nrows=1).iloc[0,0].strip().upper() == 'NULL':
        # Load descriptions into DataFrame
        df_descriptions = pd.read_csv(descriptions_path, sep='\t', header=None)
        df_descriptions.columns = ['target_id'] + [f"{db_name}_{i}" for i in range(1, len(df_descriptions.columns))]

        # Merge MMseqs output with descriptions
        df_merged = pd.merge(df_mmseqs_final, df_descriptions, left_on=f"{db_name}_id", right_on='target_id', how='left').drop('target_id', axis=1)
        df_merged.to_csv(output_path, index=False)
    else:
        # If descriptions are 'NULL', just save the processed MMseqs output
        df_mmseqs_final.to_csv(output_path, index=False)

if __name__ == "__main__":
    sample, db_name, descriptions_path, bit_score_threshold, cpus = sys.argv[1:]
    main(sample, db_name, descriptions_path, bit_score_threshold, cpus)
