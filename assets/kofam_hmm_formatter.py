import pandas as pd
import argparse
import re

def extract_ec_numbers(definition):
    """
    Extract EC numbers from the definition string and format them into a semi-colon separated list.
    Each EC number is prefixed with 'EC:'.
    """
    # This regex pattern looks for EC numbers within square brackets and captures the numbers following "EC:"
    ec_numbers = re.findall(r'\[EC:(.*?)\]', definition)
    # Join the EC numbers with semi-colon separator and prepend 'EC:' to each number
    formatted_ec_numbers = '; '.join([f"EC:{ec.strip()}" for ec_block in ec_numbers for ec in ec_block.split()])
    return formatted_ec_numbers

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
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--ch_kofam_ko", type=str, help="Path to the ch_kofam_ko file.")
    parser.add_argument("--gene_locs", type=str, help="Path to the gene locations TSV file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")
    args = parser.parse_args()

    # Load HMM search results CSV file
    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)

    # Calculate strandedness
    print("Calculating strandedness...")
    hits_df['strandedness'] = hits_df['strandedness'].apply(calculate_strandedness)

    # Load gene locations TSV file
    print("Loading gene locations TSV file...")
    gene_locs_df = pd.read_csv(args.gene_locs, sep='\t', header=None, names=['query_id', 'start_position', 'stop_position'])

    # Merge hits_df with gene_locs_df
    hits_df = pd.merge(hits_df, gene_locs_df, on='query_id', how='left')

    # Preprocess HMM search results
    print("Processing HMM search results...")
    hits_df['target_id'] = hits_df['target_id'].str.replace(r'.hmm', '', regex=True)
    hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)
    hits_df['score_rank'] = hits_df.apply(calculate_rank, axis=1)
    hits_df.dropna(subset=['score_rank'], inplace=True)

    # Find the best hit for each unique query_id
    best_hits = hits_df.groupby('query_id').apply(find_best_kofam_hit).reset_index(drop=True)

    # Mark the best hit for each unique query_id based on score_rank
    best_hits = best_hits.groupby('query_id').apply(mark_best_hit_based_on_rank).reset_index(drop=True)

    # Load ch_kofam_ko file
    print("Loading ch_kofam_ko file...")
    ch_kofam_ko_df = pd.read_csv(args.ch_kofam_ko, sep="\t")

    # Extract and format EC numbers from the "definition" column
    ch_kofam_ko_df['kofam_EC'] = ch_kofam_ko_df['definition'].apply(extract_ec_numbers)

    # Example of merging (assuming the rest of your script runs before this):
    merged_df = pd.merge(best_hits, ch_kofam_ko_df[['knum', 'definition', 'kofam_EC']], left_on='target_id', right_on='knum', how='left')

    # Keep only the relevant columns in the final output, including 'kofam_EC'
    final_output_df = merged_df[['query_id', 'start_position', 'stop_position', 'strandedness', 'target_id', 'score_rank', 'bitScore', 'definition', 'kofam_EC']]

    # Rename the columns for clarity
    final_output_df.columns = ['query_id', 'start_position', 'stop_position', 'strandedness', 'kofam_id', 'kofam_score_rank', 'kofam_bitScore', 'kofam_description', 'kofam_EC']

    # Save the modified DataFrame to CSV
    final_output_df.to_csv(args.output, index=False)

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
