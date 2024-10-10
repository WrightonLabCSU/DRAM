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
    return row['strandedness']

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
    parser.add_argument("--hits_csv", type=str, required=True, help="Path to the HMM search results CSV file.")
    parser.add_argument("--ch_canthyd_ko", type=str, required=True, help="Path to the ch_canthyd_ko file containing descriptions.")
    parser.add_argument("--gene_locs", type=str, required=True, help="Path to the gene locations TSV file.")
    parser.add_argument("--output", type=str, required=True, help="Path to the formatted output file.")
    args = parser.parse_args()

    hits_df = pd.read_csv(args.hits_csv)
    gene_locs_df = pd.read_csv(args.gene_locs, sep='\t', header=None, names=['query_id', 'start_position', 'stop_position'])
    ch_canthyd_ko_df = pd.read_csv(args.ch_canthyd_ko, sep="\t")

    # Rename the first column of ch_canthyd_ko_df to 'target_id' if it's named differently
    ch_canthyd_ko_df.rename(columns={ch_canthyd_ko_df.columns[0]: 'target_id'}, inplace=True)

    hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)
    hits_df['score_rank'] = hits_df.apply(calculate_rank, axis=1)
    hits_df['strandedness'] = hits_df.apply(calculate_strandedness, axis=1)

    # Merge hits_df with ch_canthyd_ko_df on 'target_id'
    merged_df = pd.merge(hits_df, ch_canthyd_ko_df, on='target_id', how='left')

    # Proceed with merging gene locations data and other operations
    merged_df = pd.merge(merged_df, gene_locs_df, on='query_id', how='left')

    merged_df['cant_hyd_rank'] = merged_df.apply(lambda row: assign_canthyd_rank(row, row['A_rank'], row['B_rank']), axis=1)
    merged_df.rename(columns={
        'score_rank': 'cant_hyd_score_rank',
        'bitScore': 'cant_hyd_bitScore',
        'target_id': 'cant_hyd_id',
    }, inplace=True)

    if 'description' in merged_df.columns:
        merged_df['cant_hyd_description'] = merged_df['description']
    else:
        merged_df['cant_hyd_description'] = 'No description available'

    final_output_df = merged_df[['query_id', 'start_position', 'stop_position', 'strandedness', 'cant_hyd_score_rank', 'cant_hyd_bitScore', 'cant_hyd_description', 'cant_hyd_rank', 'cant_hyd_id']]
    final_output_df.to_csv(args.output, index=False)

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
