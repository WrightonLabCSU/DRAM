import pandas as pd
import argparse
import re
from functools import partial

def get_sig_row(row, evalue_lim):
    return row['full_evalue'] < evalue_lim

def find_best_dbcan_hit(group):
    best_hit = group[group['full_evalue'] == group['full_evalue'].min()]
    return best_hit['target_id'].values[0]

def dbcan_hmmscan_formatter(hits, ch_dbcan_fam, ch_dbcan_subfam):
    # Sort hits within each group based on 'full_evalue'
    hits['rank'] = hits.groupby('query_id')['full_evalue'].rank()

    # Calculate 'score_rank'
    hits['score_rank'] = hits.groupby('query_id')['full_score'].rank(ascending=False)

    # Calculate 'bitScore'
    hits['bitScore'] = hits.groupby('query_id')['full_score'].transform('min')

    # Extract family and subfamily from target_id
    hits['family'] = hits['target_id'].str.extract(r'([A-Za-z0-9_]+)_\d*')
    hits['subfamily'] = hits['target_id'].str.extract(r'([A-Za-z0-9_]+)')

    # Extract 'dbcan-best-hit'
    hits['dbcan-best-hit'] = hits.groupby('query_id')['target_id'].transform(lambda x: x.str[:-4].unique().min())

    # Attempt to load fam_mapping file with error handling
    try:
        fam_mapping = pd.read_csv(ch_dbcan_fam, sep='\t', index_col=0, comment='#', header=None, names=['family-activities'])
    except pd.errors.EmptyDataError:
        # Handle the case when the file is empty
        print(f"The file {ch_dbcan_fam} is empty.")
        fam_mapping = pd.DataFrame(columns=['family-activities'])

    # Debug prints to check the columns in fam_mapping DataFrame
    print("Columns in fam_mapping DataFrame:", fam_mapping.columns)

    # Join 'family-activities' based on 'family'
    hits = hits.join(fam_mapping, on='family')

    # Attempt to load subfam_mapping file with error handling
    try:
        subfam_mapping = pd.read_csv(ch_dbcan_subfam, sep='\t', header=None, names=['subfam-EC', 'subfam-GenBank'])
    except pd.errors.EmptyDataError:
        # Handle the case when the file is empty
        print(f"The file {ch_dbcan_subfam} is empty.")
        subfam_mapping = pd.DataFrame(columns=['subfam-EC', 'subfam-GenBank'])

    # Debug prints to check the columns in subfam_mapping DataFrame
    print("Columns in subfam_mapping DataFrame:", subfam_mapping.columns)

    # Extract 'subfam-EC' and 'subfam-GenBank' based on 'subfamily'
    hits = hits.join(subfam_mapping.set_index(0), on='subfamily')

    # Concatenate values when there are multiple matches
    hits['subfam-EC'] = hits.groupby('query_id')['subfam-EC'].transform(lambda x: "; ".join(x))
    hits['subfam-GenBank'] = hits.groupby('query_id')['subfam-GenBank'].transform(lambda x: "; ".join(x))

    return hits[['query_id', 'target_id', 'score_rank', 'bitScore', 'family-activities', 'subfam-EC', 'subfam-GenBank', 'dbcan-best-hit']]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Format DBCAN HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--fam", type=str, help="Path to the family activities file.")
    parser.add_argument("--subfam", type=str, help="Path to the subfamily EC file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")

    args = parser.parse_args()

    # Load hits CSV
    hits_df = pd.read_csv(args.hits_csv)

    # Format hits
    formatted_hits = dbcan_hmmscan_formatter(hits_df, args.fam, args.subfam)

    # Save formatted output
    formatted_hits.to_csv(args.output, index=False)
