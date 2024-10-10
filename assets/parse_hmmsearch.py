import sys
import pandas as pd

HMMSCAN_ALL_COLUMNS = ['query_id', 'query_ascession', 'query_length', 'target_id', 'target_ascession', 'target_length', 'full_evalue', 'full_score', 'full_bias', 'domain_number', 'domain_count', 'domain_cevalue', 'domain_ievalue', 'domain_score', 'domain_bias', 'target_start', 'target_end', 'alignment_start', 'alignment_end', 'query_start', 'query_end', 'accuracy', 'description']
HMMSCAN_COLUMN_TYPES = [str, str, int, str, str, int, float, float, float, int, int, float, float, float, float, int, int, int, int, int, int, float, str]

def parse_hmmsearch_domtblout(file):
    df_lines = []
    for line in open(file):
        if not line.startswith('#'):
            line = line.split()
            line = line[:22] + [' '.join(line[22:])]
            df_lines.append(line)

    column_data = {col: [] for col in HMMSCAN_ALL_COLUMNS}
    for line in df_lines:
        for i, col in enumerate(HMMSCAN_ALL_COLUMNS):
            column_data[col].append(HMMSCAN_COLUMN_TYPES[i](line[i]))

    hmmsearch_frame = pd.DataFrame(column_data)
    
    # Extract the sign part from the description column to determine strandedness
    hmmsearch_frame['strandedness'] = hmmsearch_frame['description'].apply(lambda x: 1 if x.endswith('+') else -1)

    return hmmsearch_frame


if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    result = parse_hmmsearch_domtblout(input_file)
    result.to_csv(output_file, index=False)
