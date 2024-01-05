import pandas as pd
import argparse

def score_rank_per_row(row):
    return row['full_score'] if row['full_score'] > row['score_rank'] else row['score_rank']

def bitScore_per_row(row):
    return row['full_score'] / row['domain_number']

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format DBCAN HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--fam", type=str, help="Path to ch_dbcan_fam file.")
    parser.add_argument("--subfam", type=str, help="Path to ch_dbcan_subfam file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")

    args = parser.parse_args()

    hits_df = pd.read_csv(args.hits_csv, sep=",")  # Assuming comma as the delimiter

    # Calculate 'score_rank' and 'bitScore' for each row
    hits_df['score_rank'] = hits_df.apply(score_rank_per_row, axis=1)
    hits_df['bitScore'] = hits_df.apply(bitScore_per_row, axis=1)

    # Debugging: Print column names and contents of hits_df after calculations
    print("Column names after calculations:", hits_df.columns)
    print("Contents of hits_df after calculations:")
    print(hits_df.head())


    # Read ch_dbcan_fam file
    ch_dbcan_fam = read_ch_dbcan_fam(args.fam)

    # Read ch_dbcan_subfam file
    ch_dbcan_subfam = pd.read_csv(args.subfam, sep="\t", header=None, names=['AA', 'GenBank', 'EC'])

    # Remove ".hmm" extension from target_id in hits_df
    hits_df['target_id'] = hits_df['target_id'].apply(remove_extension)

    # Merge hits_df with ch_dbcan_fam to get subfamily
    hits_df = pd.merge(hits_df, ch_dbcan_fam, left_on='target_id', right_on='AA', how='left')

    # Merge hits_df with ch_dbcan_subfam to get GenBank and EC
    merged_df = pd.merge(hits_df, ch_dbcan_subfam, on='AA', how='left')

    # Group by query_id and aggregate multiple GenBank and EC values
    grouped_df = merged_df.groupby('query_id').agg({
        'target_id': 'first',
        'bitScore': 'first',
        'subfamily': 'first',
        'GenBank': lambda x: "; ".join(x.dropna()),  # Concatenate GenBank values
        'EC': lambda x: "; ".join(x.dropna())  # Concatenate EC values
    }).reset_index()

    # Save the formatted hits to the output file
    grouped_df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
