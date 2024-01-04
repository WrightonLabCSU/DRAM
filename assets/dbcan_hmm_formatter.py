import pandas as pd
import argparse

def dbcan_hmmscan_formater(hits, ch_dbcan_fam, ch_dbcan_subfam):
    try:
        # Debugging statement 1
        print("\nOriginal Hits DataFrame:")
        print(hits.head())

        # Check if 'query_id' column is present
        if 'query_id' not in hits.columns:
            raise KeyError("'query_id' column not found in the CSV file.")

        # Sort hits within each group based on 'full_evalue'
        hits['rank'] = hits.groupby('query_id')['full_evalue'].rank()

        # Calculate 'score_rank'
        hits['score_rank'] = hits.groupby('query_id')['full_score'].rank(ascending=False)

        # Calculate 'bitScore'
        hits['bitScore'] = hits.groupby('query_id')['full_score'].transform('min')

        # Debugging statement 2
        print("\nHits DataFrame after adding 'bitScore' column:")
        print(hits.head())

        # Extract family and subfamily from target_id
        hits['family'] = hits['target_id'].str.extract(r'([A-Za-z0-9_]+)_\d*')
        hits['subfamily'] = hits['target_id'].str.extract(r'([A-Za-z0-9_]+)')

        # Extract 'dbcan-best-hit'
        hits['dbcan-best-hit'] = hits.groupby('query_id')['target_id'].transform(lambda x: x.str[:-4].unique().min())

        # Extract 'family-activities' based on 'target_id'
        fam_mapping = pd.read_csv(ch_dbcan_fam, sep='\t', index_col=0, comment='#', header=None, names=['family-activities'], error_bad_lines=False, warn_bad_lines=False)
        fam_mapping.index = fam_mapping.index.map(str)  # Convert index to string to avoid issues
        hits['family'] = hits['family'].astype(str)  # Convert 'family' column to string
        hits = hits.join(fam_mapping, on='family')

        # Debugging statement 3
        print("\nHits DataFrame after joining with 'fam_mapping':")
        print(hits.head())

        # Extract 'subfam-EC' and 'subfam-GenBank' based on 'target_id'
        subfam_mapping = pd.read_csv(ch_dbcan_subfam, sep='\t', header=None, names=['subfam-GenBank', 'subfam-EC'])
        hits = hits.join(subfam_mapping.set_index(0), on='subfamily')

        # Debugging statement 4
        print("\nHits DataFrame after joining with 'subfam_mapping':")
        print(hits.head())

        # Concatenate values when there are multiple matches
        hits['subfam-EC'] = hits.groupby('query_id')['subfam-EC'].transform(lambda x: "; ".join(x))
        hits['subfam-GenBank'] = hits.groupby('query_id')['subfam-GenBank'].transform(lambda x: "; ".join(x))

        # Handle lines with additional columns
        additional_columns = [col for col in hits.columns[23:] if hits[col].notna().any()]
        hits['family-activities'] = hits.apply(lambda row: "; ".join(filter(None, [row['family-activities']] + row[additional_columns].astype(str))), axis=1)

        # Debugging statement 5
        print("\nHits DataFrame after final processing:")
        print(hits.head())

        # Ensure required columns are present
        required_columns = ['query_id', 'target_id', 'score_rank', 'bitScore', 'family-activities', 'subfam-EC', 'subfam-GenBank', 'dbcan-best-hit']
        selected_columns = [col for col in required_columns if col in hits.columns]

        return hits[selected_columns]

    except KeyError as e:
        print(f"Error: {e}")
        return pd.DataFrame()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Format DBCAN HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--fam", type=str, help="Path to the family activities file.")
    parser.add_argument("--subfam", type=str, help="Path to the subfamily EC file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")

    args = parser.parse_args()

    # Read CSV file with comma as the delimiter
    hits_df = pd.read_csv(args.hits_csv, delimiter=',', error_bad_lines=False)

    # Format hits
    formatted_hits = dbcan_hmmscan_formater(hits_df, args.fam, args.subfam)

    if not formatted_hits.empty:
        formatted_hits.to_csv(args.output, index=False)
        print(f"Formatted hits saved to {args.output}")
    else:
        print("Error: Unable to format hits. Please check the input data.")
