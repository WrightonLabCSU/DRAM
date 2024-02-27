import pandas as pd
import argparse

def calculate_bit_score(row):
    """Calculate bit score for each row."""
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    """Calculate rank for each row."""
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def find_best_sulfur_hit(df):
    """Find the best hit for each unique query_id."""
    # This function needs to be implemented
    return df.iloc[0]

def mark_best_hit_based_on_rank(df):
    """Mark the best hit for each unique query_id based on score_rank."""
    # This function needs to be implemented
    best_hit_idx = df["score_rank"].idxmin()
    df.at[best_hit_idx, "best_hit"] = True
    return df

def main():
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")
    args = parser.parse_args()

    # Load HMM search results CSV file
    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)

    # Preprocess HMM search results
    print("Processing HMM search results...")
    hits_df['target_id'] = hits_df['target_id'].str.replace(r'.hmm', '', regex=True)
    hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)
    hits_df['score_rank'] = hits_df.apply(calculate_rank, axis=1)
    hits_df.dropna(subset=['score_rank'], inplace=True)

    # Find the best hit for each unique query_id
    best_hits = hits_df.groupby('query_id').apply(find_best_sulfur_hit).reset_index(drop=True)

    # Mark the best hit for each unique query_id based on score_rank
    best_hits = best_hits.groupby('query_id').apply(mark_best_hit_based_on_rank).reset_index(drop=True)

    # Define sulfur_id, sulfur_score_rank, and sulfur_bitScore
    best_hits['sulfur_id'] = best_hits['target_id']
    best_hits['sulfur_score_rank'] = best_hits['score_rank']
    best_hits['sulfur_bitScore'] = best_hits['bitScore']

    # Keep only the relevant columns in the final output
    final_output_df = best_hits[['query_id', 'start_position', 'end_position', 'strandedness', 'sulfur_id', 'sulfur_score_rank', 'sulfur_bitScore', 'sulfur_definition']]

    # Save the modified DataFrame to CSV
    final_output_df.to_csv(args.output, index=False)

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
