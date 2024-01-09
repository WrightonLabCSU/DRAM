import pandas as pd
import argparse
import re

def calculate_bit_score(row):
    """Calculate bit score for each row."""
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    """Calculate rank for each row."""
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def calculate_coverage(row):
    """Calculate coverage for each row."""
    return (row['target_end'] - row['target_start']) / row['target_length']

def find_best_dbcan_hit(df):
    """Find the best hit based on E-value and coverage."""
    df['perc_cov'] = df.apply(
        lambda x: (x['target_end'] - x['target_start']) / x['target_length'], axis=1
    )
    df.sort_values(by=['full_evalue', 'perc_cov'], inplace=True, ascending=[True, False])
    return df.iloc[0]

def mark_best_hit_based_on_rank(df):
    """Mark the best hit for each unique query_id based on score_rank."""
    best_hit_idx = df["score_rank"].idxmin()
    df.at[best_hit_idx, "best_hit"] = True
    return df

def clean_ec_numbers(ec_entry):
    """Clean up EC numbers by removing '[EC:' and ']'. Replace spaces between EC numbers with ';'.

    Args:
        ec_entry (str): The input string containing EC numbers.

    Returns:
        str: The cleaned EC numbers.
    """
    ec_matches = re.findall(r'\[EC:([^\]]*?)\]', ec_entry)
    cleaned_ec_numbers = []

    for match in ec_matches:
        ec_numbers = match.split()
        for ec in ec_numbers:
            cleaned_ec = re.sub(r'[^0-9.-]', '', ec)
            cleaned_ec_numbers.append(cleaned_ec)

    result = '; '.join(cleaned_ec_numbers)
    return result

def main():
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--ch_kofam_ko", type=str, help="Path to the ch_kofam_ko file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")
    args = parser.parse_args()

    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)

    print("Processing HMM search results...")
    hits_df['target_id'] = hits_df['target_id'].str.replace(r'.hmm', '', regex=True)
    hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)
    hits_df['score_rank'] = hits_df.apply(calculate_rank, axis=1)
    hits_df['coverage'] = hits_df.apply(calculate_coverage, axis=1)
    hits_df.dropna(subset=['score_rank'], inplace=True)

    # Filter based on E-value
    hits_df = hits_df[hits_df['full_evalue'] >= 1e-18]

    # Find the best hit for each unique query_id
    best_hits = hits_df.groupby('query_id').apply(find_best_dbcan_hit).reset_index(drop=True)
    
    # Merge the best hits back to the original DataFrame
    hits_df = pd.merge(hits_df, best_hits, on='query_id', how='left')

    # Mark the best hit for each unique query_id based on score_rank
    hits_df = hits_df.groupby('query_id').apply(mark_best_hit_based_on_rank).reset_index(drop=True)

    print("Loading ch_kofam_ko file...")
    ch_kofam_ko_df = pd.read_csv(args.ch_kofam_ko, sep="\t")

    merged_df = pd.merge(hits_df, ch_kofam_ko_df[['knum', 'definition']], left_on='target_id', right_on='knum', how='left')

    merged_df['kofam_definition'] = merged_df['definition'].apply(lambda x: re.sub(r' \[EC:[^\]]*\]', '', str(x)) if pd.notna(x) else '')
    merged_df['kofam_EC'] = merged_df['definition'].apply(lambda x: clean_ec_numbers(str(x)) if pd.notna(x) else '')

    selected_columns = ['query_id', 'target_id', 'score_rank', 'bitScore', 'kofam_definition', 'kofam_EC']
    modified_columns = ['query_id', 'kofam_id', 'kofam_score_rank', 'kofam_bitScore', 'kofam_definition', 'kofam_EC']

    if set(selected_columns).issubset(merged_df.columns):
        merged_df.rename(columns=dict(zip(selected_columns, modified_columns)), inplace=True)
        merged_df[modified_columns].to_csv(args.output, index=False)

        print("Process completed successfully!")
    else:
        print("Error: Some columns are missing in the DataFrame.")

if __name__ == "__main__":
    main()
