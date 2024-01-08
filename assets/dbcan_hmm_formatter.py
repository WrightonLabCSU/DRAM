import pandas as pd
import argparse

def get_sig_row(row):
    return row['full_evalue'] < 1e-5

def calculate_bit_score(row):
    return row['full_score'] / row['domain_number']

def calculate_rank(row):
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def extract_family(target_id, ch_dbcan_fam):
    target_id_without_extension = target_id.replace('.hmm', '')
    target_id_without_underscore = target_id.split('_')[0]

    matching_rows_fam_exact = ch_dbcan_fam[ch_dbcan_fam.iloc[:, 0] == target_id_without_extension]
    matching_rows_fam_partial = ch_dbcan_fam[ch_dbcan_fam.iloc[:, 0] == target_id_without_underscore]

    matching_rows = pd.concat([matching_rows_fam_exact, matching_rows_fam_partial])

    if not matching_rows.empty:
        family_values = '; '.join(set(matching_rows.iloc[:, 1].astype(str)))  # Assuming 0-based index for columns
        family_values = family_values.strip()  # Remove leading and trailing spaces
        return family_values if pd.notna(family_values) else ""

    return ""

def extract_subfam_genbank(target_id, ch_dbcan_subfam):
    matching_rows = ch_dbcan_subfam[ch_dbcan_subfam.iloc[:, 0] == target_id]
    
    if not matching_rows.empty:
        genbank_values = '; '.join(set(matching_rows.iloc[:, 1].astype(str)))  # Assuming 0-based index for columns
        return genbank_values if pd.notna(genbank_values) else ""

    return ""

def extract_subfam_ec(target_id, ch_dbcan_subfam):
    matching_rows = ch_dbcan_subfam[ch_dbcan_subfam.iloc[:, 0] == target_id]
    
    if not matching_rows.empty:
        ec_values = '; '.join(set(matching_rows.iloc[:, 2].astype(str)))  # Assuming 0-based index for columns
        return ec_values if pd.notna(ec_values) else ""

    return ""

def find_best_dbcan_hit(df):
    df.sort_values("full_evalue", inplace=True)
    return df.iloc[0]["target_id"]


def mark_best_hit_based_on_rank(df):
    best_hit_idx = df["score_rank"].idxmin()
    df.at[best_hit_idx, "best_hit"] = True
    return df

def main():
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--fam", type=str, help="Path to the fam file.")
    parser.add_argument("--subfam", type=str, help="Path to the subfam file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")

    args = parser.parse_args()

    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)
    print("Loading subfam file...")
    ch_dbcan_subfam = pd.read_csv(args.subfam, sep="\t", comment='#', header=None, engine='python')
    print("Loading fam file...")
    ch_dbcan_fam = pd.read_csv(args.fam, comment='#', header=None, engine='python', on_bad_lines='skip', delimiter='\t', usecols=[0, 1], quoting=3)

    print("Processing HMM search results...")
    hits_df['target_id'] = hits_df['target_id'].str.replace(r'.hmm', '', regex=True)

    hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)
    hits_df['score_rank'] = hits_df.apply(calculate_rank, axis=1)
    hits_df.dropna(subset=['score_rank'], inplace=True)

    # Apply the extraction functions to create new columns
    hits_df['family'] = hits_df['target_id'].apply(lambda x: extract_family(x, ch_dbcan_fam))
    hits_df['subfam-GenBank'] = hits_df['target_id'].apply(lambda x: extract_subfam_genbank(x, ch_dbcan_subfam))
    hits_df['subfam-EC'] = hits_df['target_id'].apply(lambda x: extract_subfam_ec(x, ch_dbcan_subfam))


    # Find the best hit for each unique query_id
    best_hits = hits_df.groupby('query_id').apply(find_best_dbcan_hit).reset_index(name='dbcan-best-hit')

    # Merge the best hits back to the original DataFrame
    hits_df = pd.merge(hits_df, best_hits, on='query_id', how='left')

    # Mark the best hit for each unique query_id based on score_rank
    hits_df = hits_df.groupby('query_id').apply(mark_best_hit_based_on_rank).reset_index(drop=True)

    print("Saving the formatted output to CSV...")
    selected_columns = ['query_id', 'target_id', 'score_rank', 'bitScore', 'family', 'subfam-GenBank', 'subfam-EC', 'dbcan-best-hit']
    modified_columns = ['query_id', 'dbcan_id', 'dbcan_score_rank', 'dbcan_bitScore', 'dbcan_family', 'subfam_GenBank', 'subfam_EC', 'dbcan_best_hit']

    # Ensure the columns exist in the DataFrame before renaming
    if set(selected_columns).issubset(hits_df.columns):
        # Rename the selected columns
        hits_df.rename(columns=dict(zip(selected_columns, modified_columns)), inplace=True)
    
        # Save the modified DataFrame to CSV
        hits_df[modified_columns].to_csv(args.output, index=False)

        print("Process completed successfully!")
    else:
        print("Error: Some columns are missing in the DataFrame.")

if __name__ == "__main__":
    main()
