import pandas as pd
import argparse

def calculate_strandedness(strandedness):
    """Calculate strandedness based on the strandedness information."""
    if strandedness == '1':
        return '1'
    elif strandedness == '-1':
        return '-1'
    else:
        return strandedness

def calculate_bit_score(row):
    """Calculate bit score for each row."""
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    """Calculate rank for each row."""
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def calculate_perc_cov(row):
    """Calculate percent coverage for each row."""
    return (row['target_end'] - row['target_start']) / row['target_length']

def find_best_kofam_hit(df):
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
    parser.add_argument("--ch_camper_list", type=str, help="Path to the ch_camper_list file.")
    parser.add_argument("--gene_locs", type=str, help="Path to the gene locations TSV file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")
    args = parser.parse_args()

    # Load HMM search results CSV file
    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)

    # Load gene locations TSV file
    print("Loading gene locations TSV file...")
    gene_locs_df = pd.read_csv(args.gene_locs, sep='\t', header=None, names=['query_id', 'start_position', 'stop_position'])

    # Merge gene locations into the hits dataframe
    hits_df = pd.merge(hits_df, gene_locs_df, on='query_id', how='left')

    # Load ch_camper_list file
    print("Loading ch_camper_list file...")
    descriptions_df = pd.read_csv(args.ch_camper_list, sep="\t")

    # Rename 'target_id' in hits_df to 'camper_id' for clarity and consistency
    hits_df.rename(columns={'target_id': 'camper_id'}, inplace=True)

    # Calculate additional fields as before
    print("Calculating additional fields...")
    hits_df['strandedness'] = hits_df.apply(calculate_strandedness, axis=1)
    hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)
    hits_df['score_rank'] = hits_df.apply(calculate_rank, axis=1)
    hits_df.dropna(subset=['score_rank'], inplace=True)

    # Merge hits_df with descriptions_df
    merged_df = pd.merge(hits_df, descriptions_df, left_on='camper_id', right_on='query_id', how='left')
    merged_df.drop(columns=['query_id'], inplace=True)  # Remove duplicate 'query_id' column after merge

    # Calculate camper_rank and clean EC numbers
    print("Calculating camper_rank and cleaning EC numbers...")
    merged_df['camper_rank'] = merged_df.apply(lambda row: assign_camper_rank(row, row['A_rank'], row['B_rank']), axis=1)
    merged_df['camper_EC'] = merged_df['definition'].apply(clean_ec_numbers)

    # Finalize output, ensuring columns are correctly named and ordered
    print("Finalizing output...")
    final_columns = ['query_id', 'start_position', 'stop_position', 'strandedness', 'camper_id', 'score_rank', 'bitScore', 'camper_rank', 'camper_EC']
    final_output_df = merged_df[final_columns]
    final_output_df.columns = ['query_id', 'start_position', 'stop_position', 'strandedness', 'camper_id', 'camper_score_rank', 'camper_bitScore', 'camper_rank', 'camper_EC']

    # Save the modified DataFrame to CSV
    final_output_df.to_csv(args.output, index=False)

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
