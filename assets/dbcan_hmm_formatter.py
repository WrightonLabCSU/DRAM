import pandas as pd
import argparse

def dbcan_hmmscan_formater(hits, ch_dbcan_fam, ch_dbcan_subfam):
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

    # Extract 'family-activities' based on 'target_id'
    fam_mapping = pd.read_csv(ch_dbcan_fam, sep='\t', index_col=0, comment='#', header=None, names=['family-activities'], error_bad_lines=False)
    hits = hits.join(fam_mapping, on='family')

    # Extract 'subfam-EC' and 'subfam-GenBank' based on 'target_id'
    subfam_mapping = pd.read_csv(ch_dbcan_subfam, sep='\t', header=None, names=['subfam-EC', 'subfam-GenBank'])
    hits = hits.join(subfam_mapping.set_index(0), on='subfamily')

    # Concatenate values when there are multiple matches
    hits['subfam-EC'] = hits.groupby('query_id')['subfam-EC'].transform(lambda x: "; ".join(x))
    hits['subfam-GenBank'] = hits.groupby('query_id')['subfam-GenBank'].transform(lambda x: "; ".join(x))

    # Handle lines with additional columns
    hits['family-activities'] = hits.apply(lambda row: "; ".join(filter(None, [row['family-activities'], row.iloc[27:].astype(str)])), axis=1)

    return hits[['query_id', 'target_id', 'score_rank', 'bitScore', 'family-activities', 'subfam-EC', 'subfam-GenBank', 'dbcan-best-hit']]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Format DBCAN HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--fam", type=str, help="Path to the family activities file.")
    parser.add_argument("--subfam", type=str, help="Path to the subfamily EC file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")

    args = parser.parse_args()

    formatted_hits = dbcan_hmmscan_formater(pd.read_csv(args.hits_csv, error_bad_lines=False, delimiter='\t'), args.fam, args.subfam)
    formatted_hits.to_csv(args.output, index=False)
