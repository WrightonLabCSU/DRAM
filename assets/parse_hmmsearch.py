#!/usr/bin/env python
import sys
import pandas as pd

HMMSCAN_ALL_COLUMNS = ['query_id', 'query_ascession', 'query_length', 'target_id', 'target_ascession', 'target_length', 'full_evalue', 'full_score', 'full_bias', 'domain_number', 'domain_count', 'domain_cevalue', 'domain_ievalue', 'domain_score', 'domain_bias', 'target_start', 'target_end', 'alignment_start', 'alignment_end', 'query_start', 'query_end', 'accuracy', 'description']
HMMSCAN_COLUMN_TYPES = [str, str, int, str, str, int, float, float, float, int, int, float, float, float, float, int, int, int, int, int, int, float, str]

def parse_hmmsearch_domtblout(file):
    df_lines = list()
    for i, line in enumerate(open(file)):
        if not line.startswith('#'):
            try:
                line = line.split(maxsplit=21)
                line = line[:22] + [' '.join(line[22:])]
                df_lines.append(line)
            except Exception as e:
                print(f"Error in line {i + 1}: {e}")
                print(f"Line content: {line}")
    hmmsearch_frame = pd.DataFrame(df_lines, columns=HMMSCAN_ALL_COLUMNS)
    for i, column in enumerate(hmmsearch_frame.columns):
        hmmsearch_frame[column] = hmmsearch_frame[column].astype(HMMSCAN_COLUMN_TYPES[i])
    return hmmsearch_frame


if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    try:
        result = parse_hmmsearch_domtblout(input_file)
        result.to_csv(output_file, index=False)
    except Exception as e:
        print(f"Error: {e}")








#!/usr/bin/env python
#import sys
#import pandas as pd
#
#HMMSCAN_ALL_COLUMNS = ['query_id', 'query_ascession', 'query_length', 'target_id', 'target_ascession', 'target_length', 'full_evalue', 'full_score', 'full_bias', 'domain_number', 'domain_count', 'domain_cevalue', 'domain_ievalue', 'domain_score', 'domain_bias', 'target_start', 'target_end', 'alignment_start', 'alignment_end', 'query_start', 'query_end', 'accuracy', 'description']
#HMMSCAN_COLUMN_TYPES = [str, str, int, str, str, int, float, float, float, int, int, float, float, float, float, int, int, int, int, int, int, float, str]
#
#def parse_hmmsearch_domtblout(file):
#    df_lines = list()
#    for line in open(file):
#        if not line.startswith('#'):
#            line = line.split()
#            line = line[:22] + [' '.join(line[22:])]
#            df_lines.append(line)
#    hmmsearch_frame = pd.DataFrame(df_lines, columns=HMMSCAN_ALL_COLUMNS)
#    for i, column in enumerate(hmmsearch_frame.columns):
#        hmmsearch_frame[column] = hmmsearch_frame[column].astype(HMMSCAN_COLUMN_TYPES[i])
#    return hmmsearch_frame
#
#if __name__ == '__main__':
#    input_file = sys.argv[1]
#    output_file = sys.argv[2]
#    result = parse_hmmsearch_domtblout(input_file)
#    result.to_csv(output_file, index=False)
