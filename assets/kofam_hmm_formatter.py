import pandas as pd
import argparse
import re

def calculate_bit_score(row):
    """Calculate bit score for each row."""
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    """Calculate rank for each row."""
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def find_best_dbcan_hit(df):
    """Find the best hit based on full_evalue."""
    df.sort_values("full_evalue", inplace=True)
    return df.iloc[0]["target_id"]

def mark_best_hit_based_on_rank(df):
    """Mark the best hit for each unique query_id based on score_rank."""
    best_hit_idx = df["score_rank"].idxmin()
    df.at[best_hit_idx, "best_hit"] = True
    return df

def clean_ec_numbers(ec_entry):
    """Clean up EC numbers by removing '[EC:' and ']'."""
    ec_numbers = re.findall(r'\[EC:(.*?)\]', ec_entry)
    cleaned_ec_numbers = '; '.join(ec_numbers)
    return cleaned_ec_numbers

def main():
    # Command-line arguments
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--ch_kofam_ko", type=str, help="Path to the ch_kofam_ko file.")
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
    best_hits = hits_df.groupby('query_id').apply(find_best_dbcan_hit).reset_index(name='dbcan-best-hit')

    # Merge the best hits back to the original DataFrame
    hits_df = pd.merge(hits_df, best_hits, on='query_id', how='left')

    # Mark the best hit for each unique query_id based on score_rank
    hits_df = hits_df.groupby('query_id').apply(mark_best_hit_based_on_rank).reset_index(drop=True)

    # Load ch_kofam_ko file
    print("Loading ch_kofam_ko file...")
    ch_kofam_ko_df = pd.read_csv(args.ch_kofam_ko, sep="\t")

    # Merge hits_df with ch_kofam_ko_df
    merged_df = pd.merge(hits_df, ch_kofam_ko_df[['knum', 'definition']], left_on='target_id', right_on='knum', how='left')

    # Extract values for kofam_definition and kofam_EC
    merged_df['kofam_definition'] = merged_df['definition'].apply(lambda x: re.sub(r' \[EC:[^\]]*\]', '', str(x)) if pd.notna(x) else '')
    merged_df['kofam_EC'] = merged_df['definition'].apply(lambda x: re.search(r'\[EC:[^\]]*\]', str(x)).group(0) if re.search(r'\[EC:[^\]]*\]', str(x)) else '')

    # Clean up kofam_EC entries by removing '[EC:' and ']'
    merged_df['kofam_EC'] = merged_df['kofam_EC'].apply(clean_ec_numbers)

    # Save the formatted output to CSV
    selected_columns = ['query_id', 'target_id', 'score_rank', 'bitScore', 'kofam_definition', 'kofam_EC']
    merged_df[selected_columns].to_csv(args.output, index=False)

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
