import pandas as pd
import argparse
import re

def calculate_strandedness(row):
    """Calculate strandedness based on the strandedness information."""
    return row['strandedness']  # Assuming 'strandedness' is a column in the DataFrame

def calculate_bit_score(row):
    """Calculate bit score for each row."""
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    """Calculate rank for each row."""
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def calculate_perc_cov(row):
    """Calculate percent coverage for each row."""
    return (row['target_end'] - row['target_start']) / row['target_length']

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
    """Clean up EC numbers by adding 'EC:' prefix and formatting as a semicolon-separated list."""
    # This pattern matches individual EC numbers within the entry
    ec_matches = re.findall(r'\b\d+\.\d+\.\d+\.\-?\d*\b', ec_entry)
    
    # Prefix each EC number with 'EC:' and join them with a semicolon and space
    formatted_ec_numbers = '; '.join(['EC:' + ec for ec in ec_matches])
    
    return formatted_ec_numbers

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
    parser = argparse.ArgumentParser(description="Format HMM search results and include gene location data.")
    parser.add_argument("--hits_csv", type=str, required=True, help="Path to the HMM search results CSV file.")
    parser.add_argument("--ch_camper_list", type=str, required=True, help="Path to the ch_camper_list file.")
    parser.add_argument("--gene_locs", type=str, required=True, help="Path to the gene locations TSV file.")
    parser.add_argument("--output", type=str, required=True, help="Path to the formatted output file.")
    args = parser.parse_args()

    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)

    print("Loading gene locations TSV file...")
    gene_locs_df = pd.read_csv(args.gene_locs, sep='\t', header=None, names=['query_id', 'start_position', 'stop_position'])

    print("Loading ch_camper_list file...")
    descriptions_df = pd.read_csv(args.ch_camper_list, sep="\t")
    descriptions_df.rename(columns={'hmm_name': 'camper_id'}, inplace=True)  # Adjust based on the new information

    # IMPORTANT: Rename 'target_id' in hits_df to 'camper_id'
    hits_df.rename(columns={'target_id': 'camper_id'}, inplace=True)

    # Merge gene locations into the hits dataframe
    merged_df = pd.merge(hits_df, gene_locs_df, on='query_id', how='left')

    # Calculate additional fields
    print("Calculating additional fields...")
    merged_df['strandedness'] = merged_df.apply(calculate_strandedness, axis=1)
    merged_df['bitScore'] = merged_df.apply(calculate_bit_score, axis=1)
    merged_df['score_rank'] = merged_df.apply(calculate_rank, axis=1)
    merged_df.dropna(subset=['score_rank'], inplace=True)

    # Merge with descriptions
    merged_df = pd.merge(merged_df, descriptions_df, on='camper_id', how='left')

    print("Calculating camper_rank and cleaning EC numbers...")
    merged_df['camper_rank'] = merged_df.apply(lambda row: assign_camper_rank(row, row.get('A_rank', pd.NA), row.get('B_rank', pd.NA)), axis=1)
    merged_df['camper_EC'] = merged_df['definition'].apply(clean_ec_numbers)  # Assuming clean_ec_numbers function is defined

    print("Finalizing output...")
    final_columns = ['query_id', 'start_position', 'stop_position', 'strandedness', 'camper_id', 'score_rank', 'bitScore', 'camper_rank', 'camper_EC']
    final_output_df = merged_df[final_columns]
    final_output_df.columns = ['query_id', 'start_position', 'stop_position', 'strandedness', 'camper_id', 'camper_score_rank', 'camper_bitScore', 'camper_rank', 'camper_EC']

    final_output_df.to_csv(args.output, index=False)

    print("Process completed successfully!")

if __name__ == "__main__":
    main()