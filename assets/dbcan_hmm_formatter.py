import pandas as pd
import argparse

def get_sig_row(row):
    return row['full_evalue'] < 1e-5

def bitScore_per_row(row):
    return row['full_score'] / row['domain_number']

def rank_per_row(row):
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def process_dbcan_hits(hits, fam, subfam):
    # Create subfamily column
    hits['subfamily'] = hits['target_id'].str.rstrip('.hmm').map(fam.set_index('AA')['Description'])
    
    # Create subfam-GenBank and subfam-EC columns
    hits['subfam-GenBank'] = hits['target_id'].str.rstrip('.hmm').map(subfam.set_index('AA1_1')['GenBank'])
    hits['subfam-EC'] = hits['target_id'].str.rstrip('.hmm').map(subfam.set_index('AA1_1')['EC'])
    
    return hits[['query_id', 'target_id', 'score_rank', 'bitScore', 'subfamily', 'subfam-GenBank', 'subfam-EC']]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format KEGG HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--fam", type=str, help="Path to the ch_dbcan_fam file.")
    parser.add_argument("--subfam", type=str, help="Path to the ch_dbcan_subfam file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")

    args = parser.parse_args()

    hits = pd.read_csv(args.hits_csv)
    ch_dbcan_fam = pd.read_csv(args.fam, sep="\t", comment='#', header=None, names=['AA', 'Description'])

    # Assuming ch_dbcan_subfam has columns 'AA1_1', 'GenBank', 'EC'
    ch_dbcan_subfam = pd.read_csv(args.subfam, sep="\t", header=None, names=['AA1_1', 'GenBank', 'EC'])

    # Perform the processing
    formatted_hits = process_dbcan_hits(hits, ch_dbcan_fam, ch_dbcan_subfam)

    # Save the formatted hits
    formatted_hits.to_csv(args.output, sep="\t", index=False)
