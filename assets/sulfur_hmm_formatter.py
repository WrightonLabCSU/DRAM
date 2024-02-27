import pandas as pd
import argparse

def calculate_bit_score(row):
    """Calculate bit score for each row."""
    return row['full_score'] / row['domain_number']

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

    # Find the best hit for each unique query_id
    best_hits = hits_df.groupby('query_id').first().reset_index()

    # Add the additional columns to the output
    best_hits['start_position'] = best_hits['query_start']
    best_hits['end_position'] = best_hits['query_end']
    best_hits['strandedness'] = best_hits['strandedness']

    # Keep only the relevant columns in the final output
    final_output_df = best_hits[['query_id', 'start_position', 'end_position', 'strandedness', 'target_id', 'bitScore']]

    # Rename the columns
    final_output_df.columns = ['query_id', 'start_position', 'end_position', 'strandedness', 'sulfur_id', 'sulfur_bitScore']

    # Save the modified DataFrame to CSV
    final_output_df.to_csv(args.output, index=False)

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
