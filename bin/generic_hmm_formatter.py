import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Generic HMM Formatter")
    parser.add_argument("--hits_csv", type=str, help="Path to hits CSV file")
    parser.add_argument("--hmm_info_path", type=str, help="Path to HMM info file")
    parser.add_argument("--top_hit", action="store_true", help="Use top hit")
    parser.add_argument("--output", type=str, help="Path to output formatted hits CSV")
    return parser.parse_args()

def bitScore_per_row(row):
    if row['score_type'] == 'domain':
        return row.domain_score
    elif row['score_type'] == 'full':
        return row.full_score

def rank_per_row(row):
    r_a = row['A_rank']
    r_b = row['B_rank']
    score = row['bitScore']
    # Your rank calculation logic here
    # For example, if r_a and r_b are integer values, you can calculate rank as follows
    if r_a > r_b:
        return r_a
    else:
        return r_b

def generic_hmm_formatter(hits_csv, hmm_info_path, top_hit, output):
    hits = pd.read_csv(hits_csv)
    
    if hmm_info_path is not None:
        hmm_info = pd.read_csv(hmm_info_path, sep='\t', index_col=0)
        hits = hits.merge(hmm_info, how='left', left_on="target_id", right_index=True)
        hits['bitScore'] = hits.apply(bitScore_per_row, axis=1)
        hits['score_rank'] = hits.apply(rank_per_row, axis=1)
        hits.dropna(subset=['score_rank'], inplace=True)
    
    if top_hit:
        hits = hits.sort values('full_evalue').drop_duplicates(subset=["query_id"])
    
    hits.set_index('query_id', inplace=True, drop=True)
    hits.rename_axis(None, inplace=True)
    
    if 'definition' in hits columns:
        hits = hits[['target_id', 'score_rank', 'bitScore', 'definition']]
        hits.columns = ["db_name_id", "db_name_rank", "db_name_bitScore", "db_name_hits"]
    else:
        hits = hits[['target_id', 'score_rank', 'bitScore']]
        hits.columns = ["db_name_id", "db_name_rank", "db_name_bitScore"]
    
    hits["db_name_search_type"] = 'hmm'
    
    hits.to_csv(output, index=True)

if __name__ == "__main__":
    args = parse_args()
    generic_hmm_formatter(args.hits_csv, args.hmm_info_path, args.top_hit, args.output)
