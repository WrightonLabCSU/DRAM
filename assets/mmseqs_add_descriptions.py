import sys
import pandas as pd

def main(sample, db_name, descriptions_path, bit_score_threshold):
    # Define file paths
    input_path = f"mmseqs_out/{sample}_mmseqs_{db_name}.tsv"
    output_path = f"mmseqs_out/{sample}_mmseqs_{db_name}_formatted.csv"

    print(f"Input path: {input_path}")
    print(f"Output path: {output_path}")

    # Load MMseqs output into DataFrame
    mmseqs_columns = ['query_id', 'target_id', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df_mmseqs = pd.read_csv(input_path, sep='\t', names=mmseqs_columns)

    print("Loading MMseqs output...")

    # Add additional columns for start and end positions, and format target_id with db_name
    df_mmseqs['start_position'] = df_mmseqs['qstart']
    df_mmseqs['end_position'] = df_mmseqs['qend']
    df_mmseqs[f"{db_name}_id"] = df_mmseqs['target_id']
    df_mmseqs[f"{db_name}_bitScore"] = df_mmseqs['bitscore']
    df_mmseqs_final = df_mmseqs[['query_id', 'start_position', 'end_position', f"{db_name}_id", f"{db_name}_bitScore"]]

    print("Processing MMseqs output...")

    # Check if db_descriptions is not 'NULL'
    if descriptions_path.upper() != 'NULL':
        print("Descriptions file is not NULL. Processing...")
        df_descriptions = pd.read_csv(descriptions_path, sep='\t')
        # Rename the columns in df_descriptions with prefix
        df_descriptions = df_descriptions.add_prefix(f"{db_name}_")
        # Merge df_mmseqs_final with df_descriptions on db_name_id
        df_merged = pd.merge(df_mmseqs_final, df_descriptions, left_on=f"{db_name}_id", right_on=f"{db_name}_id", how='left')
        # Save the merged DataFrame to the output path
        df_merged.to_csv(output_path, index=False)
    else:
        print("Descriptions file is NULL. Skipping processing.")
        # Save the formatted MMseqs output to the output path
        df_mmseqs_final.to_csv(output_path, index=False)

    # Further processing if needed...

if __name__ == "__main__":
    sample, db_name, descriptions_path, bit_score_threshold = sys.argv[1:]
    main(sample, db_name, descriptions_path, bit_score_threshold)
