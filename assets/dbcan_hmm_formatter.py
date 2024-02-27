import pandas as pd
import argparse

def get_sig_row(row):
    return row['full_evalue'] < 1e-18

def calculate_bit_score(row):
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def find_best_dbcan_hit(df):
    df.sort_values(["full_evalue", "perc_cov"], inplace=True, ascending=[True, False])
    return df.iloc[0]

def mark_best_hit_based_on_rank(df):
    best_hit_idx = df["score_rank"].idxmin()
    df.at[best_hit_idx, "best_hit"] = True
    return df

def main():
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")
    args = parser.parse_args()

    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)

    print("Processing HMM search results...")
    hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)
    hits_df['score_rank'] = hits_df.apply(calculate_rank, axis=1)
    hits_df['perc_cov'] = (hits_df['target_end'] - hits_df['target_start']) / hits_df['target_length']
    hits_df.dropna(subset=['score_rank'], inplace=True)

    # Extract start_position, end_position, and strandedness
    hits_df['start_position'] = hits_df['query_start']
    hits_df['end_position'] = hits_df['query_end']
    hits_df['strandedness'] = hits_df['strandedness']

    # Filter based on E-value
    hits_df = hits_df[hits_df.apply(get_sig_row, axis=1)]

    # Find the best hit for each unique query_id
    best_hits = hits_df.groupby('query_id').apply(find_best_dbcan_hit)

    # Keep only the rows with the best hits
    hits_df = best_hits.reset_index(drop=True)

    # Mark the best hit for each unique query_id based on score_rank
    hits_df = hits_df.groupby('query_id').apply(mark_best_hit_based_on_rank).reset_index(drop=True)

    print("Saving the formatted output to CSV...")
    selected_columns = ['query_id', 'start_position', 'end_position', 'strandedness', 'target_id', 'score_rank', 'bitScore']
    modified_columns = ['query_id', 'start_position', 'end_position', 'strandedness', 'dbcan_id', 'dbcan_score_rank', 'dbcan_bitScore']

    # Ensure the columns exist in the DataFrame before renaming
    if set(selected_columns).issubset(hits_df.columns):
        # Rename the selected columns
        hits_df.rename(columns=dict(zip(selected_columns, modified_columns)), inplace=True)

        # Save the formatted output to CSV
        hits_df[modified_columns].to_csv(args.output, index=False)

        print("Process completed successfully!")

if __name__ == "__main__":
    main()
