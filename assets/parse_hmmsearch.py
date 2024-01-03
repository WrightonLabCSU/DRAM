#!/usr/bin/env python
import sys
import pandas as pd

HMMSCAN_ALL_COLUMNS = ['query_id', 'query_ascession', 'query_length', 'target_id', 'target_ascession', 'target_length', 'full_evalue', 'full_score', 'full_bias', 'domain_number', 'domain_count', 'domain_cevalue', 'domain_ievalue', 'domain_score', 'domain_bias', 'target_start', 'target_end', 'alignment_start', 'alignment_end', 'query_start', 'query_end', 'accuracy', 'description']
HMMSCAN_COLUMN_TYPES = [str, str, int, str, str, int, float, float, float, int, int, float, float, float, float, int, int, int, int, int, int, float, str]

def parse_hmmsearch_domtblout(file):
    df_lines = []
    for line in open(file):
        if not line.startswith('#'):
            line = line.split(maxsplit=21)
            line = line[:22] + [' '.join(line[22:])]
            df_lines.append(line)

    if not df_lines:
        print("No valid lines found in the input file.")
        return pd.DataFrame()

    try:
        hmmsearch_frame = pd.DataFrame(df_lines, columns=HMMSCAN_ALL_COLUMNS)
        for i, column in enumerate(hmmsearch_frame.columns):
            hmmsearch_frame[column] = hmmsearch_frame[column].astype(HMMSCAN_COLUMN_TYPES[i])
        return hmmsearch_frame
    except Exception as e:
        print(f"Error in creating DataFrame: {e}")
        print("DataFrame creation failed. Lines:")
        for line in df_lines:
            print(line)
        return pd.DataFrame()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    try:
        result = parse_hmmsearch_domtblout(input_file)
        result.to_csv(output_file, index=False)
    except Exception as e:
        print(f"Error: {e}")
