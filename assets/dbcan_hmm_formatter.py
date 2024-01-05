import pandas as pd

def generate_subfam_genbank(row, ch_dbcan_subfam):
    matching_rows = ch_dbcan_subfam[ch_dbcan_subfam['subfamily'] == row['subfam-CAZy']]
    selected_row = matching_rows.loc[matching_rows['score'].idxmax()]
    return f"{row['subfam-CAZy']}_{selected_row['GenBank']}"

def main():
    # Replace 'path_to_ch_dbcan_fam' and 'path_to_ch_dbcan_subfam' with the actual paths to your files
    ch_dbcan_fam = pd.read_csv('path_to_ch_dbcan_fam', comment='#', header=None, names=['target_id', 'subfamily'], engine='python', error_bad_lines=False, delimiter='\t', usecols=[0, 1], quoting=3, na_values=['nan'])
    ch_dbcan_subfam = pd.read_csv('path_to_ch_dbcan_subfam', delimiter='\t')

    # Check for NaN values in the 'target_id' column of ch_dbcan_fam
    nan_values = ch_dbcan_fam['target_id'].isna().sum()
    
    if nan_values > 0:
        print(f"Warning: Found {nan_values} NaN values in 'target_id' column of ch_dbcan_fam. Consider handling NaN values.")

    # Assuming you want to exclude rows with NaN values, you can drop them
    ch_dbcan_fam = ch_dbcan_fam.dropna(subset=['target_id'])

    # Assuming you want to replace NaN values in 'target_id' with a default value (e.g., 'unknown_target')
    # ch_dbcan_fam['target_id'].fillna('unknown_target', inplace=True)

    hits_df = pd.DataFrame(...)  # Replace ... with your code to create or load hits_df

    hits_df['subfam-GenBank'] = hits_df.apply(lambda row: generate_subfam_genbank(row, ch_dbcan_subfam), axis=1)

    hits_df['subfam-EC'] = hits_df.apply(lambda row: generate_subfam_ec(row, ch_dbcan_subfam), axis=1)

    # Filter significant rows
    sig_hits_df = hits_df[hits_df.apply(get_sig_row, axis=1)]

    # Sort the significant hits by score rank
    sig_hits_df = sig_hits_df.sort_values(by='score_rank')

    # Save the formatted output to a new CSV file
    selected_columns = ['query_id', 'target_id', 'score_rank', 'bitScore', 'subfamily', 'subfam-GenBank', 'subfam-EC']
    sig_hits_df[selected_columns].to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
