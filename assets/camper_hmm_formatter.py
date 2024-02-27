import pandas as pd
import argparse
import re

def calculate_bit_score(row):
    """Calculate bit score for each row."""
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    """Calculate rank for each row."""
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def find_best_camper_hit(df):
    """Find the best hit based on E-value and coverage."""
    df['perc_cov'] = (df['target_end'] - df['target_start']) / df['target_length']
    df.sort_values(by=["full_evalue", "perc_cov"], ascending=[True, False], inplace=True)
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
    cleaned_ec_numbers = [re.sub(r'[^0-9.-]', '', ec) for match in ec_matches for ec in match.split()]
    result = '; '.join(cleaned_ec_numbers)
    return result

def assign_camper_rank(row, a_rank, b_rank):
    """Assign camper rank based on bit score and provided thresholds."""
    if pd.isna(row['bitScore']):
        return None
    elif row['bitScore'] >= a_rank:
        return 'A'
    elif row['bitScore'] >= b_rank:
        return 'B'
    else:
        return None




def main():
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--ch_camper_list", type=str, help="Path to the ch_camper_list file.")
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
    best_hits = hits_df.groupby('query_id').apply(find_best_camper_hit).reset_index(drop=True)

    # Mark the best hit for each unique query_id based on score_rank
    best_hits = best_hits.groupby('query_id').apply(mark_best_hit_based_on_rank).reset_index(drop=True)

    # Load ch_camper_list file
    print("Loading ch_camper_list file...")
    ch_camper_list_df = pd.read_csv(args.ch_camper_list, sep="\t")

    # Merge hits_df with ch_camper_list_df
    merged_df = pd.merge(best_hits, ch_camper_list_df[['hmm_name', 'A_rank', 'B_rank', 'score_type', 'definition']], left_on='target_id', right_on='hmm_name', how='left')

    # Extract values for camper_definition
    merged_df['camper_definition'] = merged_df['definition']

    # Calculate camper_rank
    merged_df['camper_rank'] = merged_df.apply(lambda row: assign_camper_rank(row, row['A_rank'], row['B_rank']), axis=1)

    # Add the additional columns to the output
    merged_df['start_position'] = merged_df['query_start']
    merged_df['end_position'] = merged_df['query_end']
    merged_df['strandedness'] = merged_df['strandedness']

    # Keep only the relevant columns in the final output
    final_output_df = merged_df[['query_id', 'start_position', 'end_position', 'strandedness', 'target_id', 'score_rank', 'bitScore', 'camper_definition', 'camper_rank']]

    # Rename the columns
    final_output_df.columns = ['query_id', 'start_position', 'end_position', 'strandedness', 'camper_id', 'camper_score_rank', 'camper_bitScore', 'camper_definition', 'camper_rank']

    # Perform EC cleaning
    final_output_df['camper_EC'] = final_output_df['camper_definition'].apply(clean_ec_numbers)

    # Save the modified DataFrame to CSV
    final_output_df.to_csv(args.output, index=False)

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
