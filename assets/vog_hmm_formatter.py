import pandas as pd
import argparse

def calculate_bit_score(row):
    """Calculate bit score for each row."""
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    """Calculate rank for each row."""
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def calculate_strandedness(row):
    """Calculate strandedness based on the strandedness information."""
    return row['strandedness']  # Assuming 'strandedness' is a column in the DataFrame


def find_best_vog_hit(df):
    """Find the best hit based on E-value and coverage."""
    df['perc_cov'] = (df['target_end'] - df['target_start']) / df['target_length']
    df.sort_values(by=["full_evalue", "perc_cov"], ascending=[True, False], inplace=True)
    return df.iloc[0]

def mark_best_hit_based_on_rank(df):
    """Mark the best hit for each unique query_id based on score_rank."""
    best_hit_idx = df["score_rank"].idxmin()
    df.at[best_hit_idx, "best_hit"] = True
    return df

def main():
    parser = argparse.ArgumentParser(description="Format HMM search results and include gene location data.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--gene_locs", type=str, help="Path to the gene locations TSV file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")
    args = parser.parse_args()

    # Load HMM search results CSV file
    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)

    # Load gene locations TSV file
    print("Loading gene locations TSV file...")
    gene_locs_df = pd.read_csv(args.gene_locs, sep='\t', header=None, names=['query_id', 'start_position', 'stop_position'])

    # Merge hits_df with gene_locs_df on query_id
    merged_df = pd.merge(hits_df, gene_locs_df, on='query_id', how='left')

    # Calculate strandedness
    print("Calculating strandedness...")
    hits_df['strandedness'] = hits_df.apply(calculate_strandedness, axis=1)

    # Process HMM search results
    print("Processing HMM search results...")
    merged_df['target_id'] = merged_df['target_id'].str.replace(r'.hmm', '', regex=True)
    merged_df['bitScore'] = merged_df.apply(calculate_bit_score, axis=1)
    merged_df['score_rank'] = merged_df.apply(calculate_rank, axis=1)
    merged_df.dropna(subset=['score_rank'], inplace=True)

    # Find and mark best hits
    best_hits = merged_df.groupby('query_id').apply(find_best_vog_hit).reset_index(drop=True)
    best_hits = best_hits.groupby('query_id').apply(mark_best_hit_based_on_rank).reset_index(drop=True)

    # Keep only the relevant columns in the final output
    final_output_df = best_hits[['query_id', 'start_position', 'stop_position', 'strandedness', 'target_id', 'score_rank', 'bitScore']]

    # Rename the columns for the final output
    final_output_df.columns = ['query_id', 'start_position', 'stop_position', 'strandedness', 'vogdb_id', 'vog_score_rank', 'vog_bitScore']

    # Save the modified DataFrame to CSV
    final_output_df.to_csv(args.output, index=False)

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
