import pandas as pd
import argparse

def calculate_bit_score(row):
    """Calculate bit score for each row."""
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    """Calculate rank for each row."""
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def assign_canthyd_rank(row, a_rank, b_rank):
    """Assign canthyd rank based on bit score and provided thresholds."""
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
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--ch_canthyd_ko", type=str, help="Path to the ch_canthyd_ko file containing the descriptions.")
    parser.add_argument("--gene_locs", type=str, help="Path to the gene locations TSV file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")
    args = parser.parse_args()

    # Load HMM search results CSV file
    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)

    # Load gene locations TSV file
    print("Loading gene locations TSV file...")
    gene_locs_df = pd.read_csv(args.gene_locs, sep='\t', header=None, names=['query_id', 'start_position', 'stop_position'])

    # Preprocess HMM search results
    print("Processing HMM search results...")
    hits_df['target_id'] = hits_df['target_id'].str.replace(r'.hmm', '', regex=True)
    hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)
    hits_df['score_rank'] = hits_df.apply(calculate_rank, axis=1)
    hits_df.dropna(subset=['score_rank'], inplace=True)

    # Load ch_canthyd_ko file
    print("Loading ch_canthyd_ko file...")
    ch_canthyd_ko_df = pd.read_csv(args.ch_canthyd_ko, sep="\t")

    # Merging HMM search results with canthyd ko data
    merged_df = pd.merge(hits_df, ch_canthyd_ko_df, left_on='target_id', right_on='hmm_name', how='left')

    # Merge gene locations data
    merged_df = pd.merge(merged_df, gene_locs_df, on='query_id', how='left')

    # Assign canthyd rank based on bit score
    merged_df['canthyd_rank'] = merged_df.apply(lambda row: assign_canthyd_rank(row, row['A_rank'], row['B_rank']), axis=1)

    # Rename 'description' column to 'canthyd_description'
    merged_df.rename(columns={'description': 'canthyd_description'}, inplace=True)

    # Prepare final output DataFrame
    final_output_df = merged_df[['query_id', 'start_position', 'stop_position', 'strandedness', 'target_id', 'score_rank', 'bitScore', 'canthyd_description', 'canthyd_rank']]

    # Rename columns appropriately
    final_output_df.columns = ['query_id', 'start_position', 'stop_position', 'strandedness', 'canthyd_id', 'canthyd_score_rank', 'canthyd_bitScore', 'canthyd_description', 'canthyd_rank']

    # Save the modified DataFrame to CSV
    final_output_df.to_csv(args.output, index=False)

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
