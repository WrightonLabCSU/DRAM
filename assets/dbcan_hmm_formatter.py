import pandas as pd
import argparse

def get_sig_row(row):
    return row['full_evalue'] < 1e-5

def bitScore_per_row(row):
    return row['full_score'] / row['domain_number']

def rank_per_row(row):
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def generate_subfamily(row, ch_dbcan_subfam):
    target_id = row['target_id']
    matching_rows = ch_dbcan_subfam[ch_dbcan_subfam['target_id'] == target_id]
    if not matching_rows.empty:
        return "; ".join(matching_rows['subfamily'])
    else:
        return ""

def generate_subfam_GenBank(row, ch_dbcan_subfam):
    target_id = row['target_id']
    matching_rows = ch_dbcan_subfam[ch_dbcan_subfam['target_id'] == target_id]
    if not matching_rows.empty:
        return "; ".join(matching_rows['subfam-GenBank'])
    else:
        return ""

def generate_subfam_EC(row, ch_dbcan_subfam):
    target_id = row['target_id']
    matching_rows = ch_dbcan_subfam[ch_dbcan_subfam['target_id'] == target_id]
    if not matching_rows.empty:
        return "; ".join(matching_rows['subfam-EC'])
    else:
        return ""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--fam", type=str, help="Path to the fam file.")
    parser.add_argument("--subfam", type=str, help="Path to the subfam file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")

    args = parser.parse_args()

    hits_df = pd.read_csv(args.hits_csv)
    ch_dbcan_subfam = pd.read_csv(args.subfam, sep="\t", comment='#', header=None, names=['target_id', 'subfamily', 'subfam-GenBank', 'subfam-EC'])

    print("Contents of ch_dbcan_subfam:")
    print(ch_dbcan_subfam.head())


    hits_df['bitScore'] = hits_df.apply(bitScore_per_row, axis=1)
    hits_df['score_rank'] = hits_df.apply(rank_per_row, axis=1)
    hits_df.dropna(subset=['score_rank'], inplace=True)

    hits_df['subfamily'] = hits_df.apply(lambda row: generate_subfamily(row, ch_dbcan_subfam), axis=1)
    hits_df['subfam-GenBank'] = hits_df.apply(lambda row: generate_subfam_GenBank(row, ch_dbcan_subfam), axis=1)
    hits_df['subfam-EC'] = hits_df.apply(lambda row: generate_subfam_EC(row, ch_dbcan_subfam), axis=1)

    print("Column names of hits_df:", hits_df.columns)
    print("Contents of hits_df:")
    print(hits_df)

    hits_df = hits_df[['query_id', 'target_id', 'score_rank', 'bitScore', 'subfamily', 'subfam-GenBank', 'subfam-EC']]
    hits_df.to_csv(args.output, sep="\t", index=False)
