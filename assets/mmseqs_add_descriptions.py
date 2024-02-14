import sys
import pandas as pd

def main(sample, db_name, descriptions_path, bit_score_threshold):
    # Define file paths
    input_path = f"mmseqs_out/{sample}_mmseqs_{db_name}.tsv"
    output_path = f"mmseqs_out/{sample}_mmseqs_{db_name}_formatted.csv"

    print("Input path:", input_path)
    print("Output path:", output_path)

    # Load MMseqs output into DataFrame
    print("Loading MMseqs output...")
    mmseqs_columns = ['query_id', 'target_id', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df_mmseqs = pd.read_csv(input_path, sep='\t', names=mmseqs_columns)

    # Add additional columns for start and end positions, and format target_id with db_name
    print("Processing MMseqs output...")
    df_mmseqs['start_position'] = df_mmseqs['qstart']
    df_mmseqs['end_position'] = df_mmseqs['qend']
    df_mmseqs[f"{db_name}_id"] = df_mmseqs['target_id']
    df_mmseqs[f"{db_name}_bitScore"] = df_mmseqs['bitscore']
    df_mmseqs_final = df_mmseqs[['query_id', 'start_position', 'end_position', f"{db_name}_id", f"{db_name}_bitScore"]]

    # Check if db_descriptions is not 'NULL'
    print("Checking descriptions...")
    with open(descriptions_path, 'r') as file:
        first_line = next(file).strip()
        print("First line of descriptions:", first_line)
        if first_line.upper() != 'NULL':
            print("Descriptions file is not NULL. Processing...")
            file.seek(0)
            df_descriptions = pd.read_csv(descriptions_path, sep='\t', header=None)
            # Process descriptions and merge...
            # Rest of the processing as before...
        else:
            print("Descriptions file is NULL. Saving MMseqs output...")
            df_mmseqs.to_csv(output_path, index=False)

if __name__ == "__main__":
    sample, db_name, descriptions_path, bit_score_threshold = sys.argv[1:]
    print("Sample:", sample)
    print("DB Name:", db_name)
    print("Descriptions Path:", descriptions_path)
    print("Bit Score Threshold:", bit_score_threshold)
    main(sample, db_name, descriptions_path, bit_score_threshold)
