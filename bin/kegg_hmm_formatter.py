import pandas as pd
import argparse

def get_sig_row(row):
    return row['full_evalue'] < 1e-5

def bitScore_per_row(row):
    return row['full_score'] / row['domain_number']

def rank_per_row(row):
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format KEGG HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--hmm_info_path", type=str, help="Path to the KEGG HMM info TSV file.")
    parser.add_argument("--top_hit", type=bool, help="Select only the top hit for each query.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")

    args = parser.parse_args()

    if args.hmm_info_path == "empty":
        # If hmm_info_path is "/empty," assume no KEGG HMM info is available.
        hits = pd.read_csv(args.hits_csv)
        hits['bitScore'] = hits.apply(bitScore_per_row, axis=1)
        hits['score_rank'] = hits.apply(rank_per_row, axis=1)
        hits.dropna(subset=['score_rank'], inplace=True)
        hits = hits[['query_id', 'target_id', 'score_rank', 'bitScore']]
    else:
        hits = pd.read_csv(args.hits_csv)
        hmm_info = pd.read_csv(args.hmm_info_path, sep="\t", index_col=0)

        if args.top_hit:
            hits = hits[hits.groupby('query_id')['full_evalue'].transform(min) == hits['full_evalue']]

        hits['bitScore'] = hits.apply(bitScore_per_row, axis=1)
        hits['score_rank'] = hits.apply(rank_per_row, axis=1)
        hits.dropna(subset=['score_rank'], inplace=True)

        hits = hits[['query_id', 'target_id', 'score_rank', 'bitScore']]
        if 'definition' in hits:
            hits = hits[['query_id', 'target_id', 'score_rank', 'bitScore', 'definition']]

    hits.to_csv(args.output, sep="\t", index=False)
