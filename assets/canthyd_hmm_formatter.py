import pandas as pd
import argparse

def calculate_bit_score(row):
    """Calculate bit score for each row."""
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    """Calculate rank for each row."""
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def main():
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--ch_canthyd_ko", type=str, help="Path to the ch_canthyd_ko file.")
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

    # Load ch_canthyd_ko file
    print("Loading ch_canthyd_ko file...")
    ch_canthyd_ko_df = pd.read_csv(args.ch_canthyd_ko, sep="\t")

    # Check if 'description' column exists in ch_canthyd_ko_df
    if 'description' not in ch_canthyd_ko_df.columns:
        print("Error: 'description' column not found in ch_canthyd_ko file.")
        return

    # Check if 'hmm_name' column exists in ch_canthyd_ko_df
    if 'hmm_name' not in ch_canthyd_ko_df.columns:
        print("Error: 'hmm_name' column not found in ch_canthyd_ko file.")
        return

    # Merge hits_df with ch_canthyd_ko_df
    print("Merging dataframes...")
    merged_df = pd.merge(hits_df, ch_canthyd_ko_df, left_on='target_id', right_on='hmm_name', how='left')

    # Extract values for canthyd_description
    merged_df['canthyd_description'] = merged_df['description']

    # Add the additional columns to the output
    merged_df['start_position'] = merged_df['query_start']
    merged_df['end_position'] = merged_df['query_end']
    merged_df['strandedness'] = merged_df['strandedness']

    # Keep only the relevant columns in the final output
    final_output_df = merged_df[['query_id', 'start_position', 'end_position', 'strandedness', 'target_id', 'score_rank', 'bitScore', 'canthyd_description']]

    # Rename the columns
    final_output_df.columns = ['query_id', 'start_position', 'end_position', 'strandedness', 'canthyd_id', 'canthyd_score_rank', 'canthyd_bitScore', 'canthyd_description']

    # Save the modified DataFrame to CSV
    final_output_df.to_csv(args.output, index=False)

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
