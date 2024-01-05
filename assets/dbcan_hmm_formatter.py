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
    return "; ".join(matching_rows['subfamily']) if not matching_rows.empty else ""

def generate_subfam_GenBank(row, ch_dbcan_subfam):
    target_id = row['target_id']
    matching_rows = ch_dbcan_subfam[ch_dbcan_subfam['target_id'] == target_id]
    return "; ".join(matching_rows['subfam-GenBank']) if not matching_rows.empty else ""

def generate_subfam_EC(row, ch_dbcan_subfam):
    target_id = row['target_id']
    matching_rows = ch_dbcan_subfam[ch_dbcan_subfam['target_id'] == target_id]
    return "; ".join(map(str, matching_rows['subfam-EC'].dropna())) if not matching_rows.empty else ""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--fam", type=str, help="Path to the fam file.")
    parser.add_argument("--subfam", type=str, help="Path to the subfam file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")

    args = parser.parse_args()

    # Read HMM search results CSV file and subfam file
    hits_df = pd.read_csv(args.hits_csv)
    ch_dbcan_subfam = pd.read_csv(args.subfam, sep="\t", comment='#', header=None, names=['target_id', 'subfamily', 'subfam-GenBank', 'subfam-EC'])

    # Remove the '.hmm' extension from 'target_id' in hits_df
    hits_df['target_id'] = hits_df['target_id'].str.replace(r'.hmm', '', regex=True)

    # Print unique target_id values in hits_df after modification
    print("\nUnique target_id values in hits_df after modification:")
    print(hits_df['target_id'].unique().tolist())

    # Print unique target_id values in ch_dbcan_subfam
    print("\nUnique target_id values in ch_dbcan_subfam:")
    print(ch_dbcan_subfam['target_id'].unique().tolist())

    # Display the target_id structure in hits_df
    print("\nStructure of target_id values in hits_df:")
    print(hits_df['target_id'].head(10))

    # Display the target_id structure in ch_dbcan_subfam
    print("\nStructure of target_id values in ch_dbcan_subfam:")
    print(ch_dbcan_subfam['target_id'].head(10))

    print("Contents of ch_dbcan_subfam:")
    print(ch_dbcan_subfam.head())

    # Add new columns to hits_df
    hits_df['bitScore'] = hits_df.apply(bitScore_per_row, axis=1)
    hits_df['score_rank'] = hits_df.apply(rank_per_row, axis=1)
    hits_df.dropna(subset=['score_rank'], inplace=True)

    # Filter matching rows between hits_df and ch_dbcan_subfam
    matching_rows = hits_df[hits_df['target_id'].isin(ch_dbcan_subfam['target_id'])]
    print("Matching rows between hits_df and ch_dbcan_subfam:")
    print(matching_rows[['query_id', 'target_id', 'score_rank', 'bitScore']])

    # Remove duplicates from ch_dbcan_subfam DataFrame
    ch_dbcan_subfam = ch_dbcan_subfam.drop_duplicates(subset='target_id')

    # Update the mapping in hits_df with correct column assignments
    hits_df['subfamily'] = hits_df['target_id'].map(ch_dbcan_subfam.set_index('target_id')['subfamily'])
    hits_df['subfam-GenBank'] = hits_df['target_id'].map(ch_dbcan_subfam.set_index('target_id')['subfam-GenBank'])
    hits_df['subfam-EC'] = hits_df['target_id'].map(ch_dbcan_subfam.set_index('target_id')['subfam-EC'])

    # Print column names and contents of hits_df
    print("Column names of hits_df:", hits_df.columns)
    print("Contents of hits_df:")
    print(hits_df)

    # Save the formatted output to a file
    hits_df[['query_id', 'target_id', 'score_rank', 'bitScore', 'subfamily', 'subfam-GenBank', 'subfam-EC']].to_csv(args.output, sep="\t", index=False)
