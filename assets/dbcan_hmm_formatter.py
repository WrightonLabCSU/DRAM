import pandas as pd
import argparse

def extract_subfamily_description(row, fam_data):
    target_id = row['target_id'].replace('.hmm', '')  # Remove the ".hmm" extension
    if target_id in fam_data.index:
        return fam_data.loc[target_id, 'Description']
    else:
        return None

def extract_genbank_ec(row, subfam_data):
    target_id = row['target_id'].replace('.hmm', '')  # Remove the ".hmm" extension
    if target_id in subfam_data.index:
        return subfam_data.loc[target_id, 'Genbank'], subfam_data.loc[target_id, 'EC']
    else:
        return None, None

def get_sig_row(row):
    return row['full_evalue'] < 1e-5

def bitScore_per_row(row):
    return row['full_score'] / row['domain_number']

def rank_per_row(row):
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format KEGG HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--fam", type=str, help="Path to the CAZyDB fam-activities file.")
    parser.add_argument("--subfam", type=str, help="Path to the CAZyDB fam.subfam.ec file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")

    args = parser.parse_args()

    fam_data = pd.read_csv(args.fam, sep="\t", comment='#', header=None, names=['AA', 'Description'], usecols=[0, 1], index_col=0)
    subfam_data = pd.read_csv(args.subfam, sep="\t", header=None, names=['Subfamily', 'Genbank', 'EC'], usecols=[0, 1, 2], index_col=0)
    hits = pd.read_csv(args.hits_csv)

    hits['subfamily'] = hits.apply(lambda row: extract_subfamily_description(row, fam_data), axis=1)
    hits['subfam-GenBank'], hits['subfam-EC'] = zip(*hits.apply(lambda row: extract_genbank_ec(row, subfam_data), axis=1))
    hits['bitScore'] = hits.apply(bitScore_per_row, axis=1)
    hits['score_rank'] = hits.apply(rank_per_row, axis=1)
    hits.dropna(subset=['score_rank'], inplace=True)
    hits = hits[['query_id', 'target_id', 'score_rank', 'bitScore', 'subfamily', 'subfam-GenBank', 'subfam-EC']]

    hits.to_csv(args.output, sep="\t", index=False)
