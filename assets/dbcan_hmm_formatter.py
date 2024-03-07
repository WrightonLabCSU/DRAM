import pandas as pd
import argparse

def get_sig_row(row):
    """Filter rows based on a significance threshold for the e-value."""
    return row['full_evalue'] < 1e-18

def calculate_bit_score(row):
    """Calculate bit score for each row."""
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    """Calculate rank for each row."""
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def find_best_dbcan_hit(df):
    """Find the best DBCAN hit based on E-value and coverage."""
    df.sort_values(["full_evalue", "perc_cov"], inplace=True, ascending=[True, False])
    return df.iloc[0]

def mark_best_hit_based_on_rank(df):
    """Mark the best hit for each unique query_id based on score_rank."""
    best_hit_idx = df["score_rank"].idxmin()
    df.at[best_hit_idx, "best_hit"] = True
    return df

def calculate_strandedness(row):
    """Calculate strandedness based on the 'strandedness' column."""
    strandedness = row['strandedness']
    if strandedness in ['+', '-']:
        return strandedness
    else:
        return ''

def main():
    parser = argparse.ArgumentParser(description="Format HMM search results and include gene location data.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--gene_locs", type=str, help="Path to the gene locations TSV file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")
    args = parser.parse_args()

    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)

    print("Loading gene locations TSV file...")
    gene_locs_df = pd.read_csv(args.gene_locs, sep='\t', header=None, names=['query_id', 'start_position', 'stop_position'])

    print("Processing HMM search results...")
    # Merge gene locations into the hits dataframe
    hits_df = pd.merge(hits_df, gene_locs_df, on='query_id', how='left')

    print("Calculating strandedness...")
    # Calculate strandedness based on the strandedness column
    hits_df['strandedness'] = hits_df.apply(calculate_strandedness, axis=1)

    # Calculate bit score
    hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)

    # Calculate rank
    hits_df['score_rank'] = hits_df.apply(calculate_rank, axis=1)

    # Filter significant rows
    hits_df = hits_df[hits_df.apply(get_sig_row, axis=1)]

    # Find the best DBCAN hit
    hits_df = find_best_dbcan_hit(hits_df)

    # Mark the best hit based on rank
    hits_df = mark_best_hit_based_on_rank(hits_df)

    # Save the formatted output to CSV including strandedness
    selected_columns = ['query_id', 'start_position', 'stop_position', 'strandedness', 'target_id', 'score_rank', 'bitScore']
    modified_columns = ['query_id', 'start_position', 'stop_position', 'strandedness', 'dbcan_id', 'dbcan_score_rank', 'dbcan_bitScore']

    # Ensure the columns exist in the DataFrame before renaming
    if set(selected_columns).issubset(hits_df.columns):
        # Rename the selected columns to include strandedness
        hits_df.rename(columns=dict(zip(selected_columns, modified_columns)), inplace=True)

        # Remove '.hmm' extension from 'dbcan_id' values if it exists
        hits_df['dbcan_id'] = hits_df['dbcan_id'].str.replace('.hmm', '', regex=False)

        # Save the formatted output to CSV including strandedness
        hits_df[modified_columns].to_csv(args.output, index=False)

        print("Process completed successfully!")

if __name__ == "__main__":
    main()
