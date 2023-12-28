import pandas as pd
import argparse
import re
from functools import partial

def get_sig_row(row, evalue_lim):
    return row['full_evalue'] < evalue_lim

def find_best_dbcan_hit(group):
    best_hit = group[group['full_evalue'] == group['full_evalue'].min()]
    return best_hit['target_id'].values[0]

def dbcan_hmmscan_formater(hits, db_name, db_handler=None):
    hits_sig = hits[hits.apply(partial(get_sig_row, evalue_lim=1e-18), axis=1]

    if len(hits_sig) == 0:
        return pd.DataFrame()

    hit_groups = hits_sig.groupby('query_id')
    all_hits = hit_groups.apply(lambda x: "; ".join(x['target_id'].apply(lambda y: y[:-4]).unique()))

    hits_df = pd.DataFrame(all_hits)
    hits_df.columns = [f"{db_name}_ids"]

    def description_pull(x):
        id_list = [re.findall("^[A-Z]*[0-9]*", str(x))[0] for x in x.split("; ")]
        id_list = [y for x in id_list for y in x if len(x) > 0]
        description_list = db_handler.get_descriptions(id_list, "dbcan_description").values()
        description_str = "; ".join(description_list)
        return description_str

    if db_handler is not None:
        hits_df[f"{db_name}_hits"] = hits_df[f"{db_name}_ids"].apply(description_pull)
        hits_df[f"{db_name}_subfam_ec"] = hits_df[f"{db_name}_ids"].apply(lambda x: "; ".join(
            db_handler.get_descriptions(x.split("; "), "dbcan_description", description_name="ec").values()))

    hits_df[f"{db_name}_best_hit"] = [find_best_dbcan_hit(group) for _, group in hit_groups]
    hits_df.rename_axis(None, inplace=True)

    return hits_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Format DBCAN HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--db_name", type=str, help="Name of the DBCAN database.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")

    args = parser.parse_args()

    # Mock database handler (replace with actual data)
    db_handler = {
        "EC1": "Description for EC1",
        "EC2": "Description for EC2",
        # Add more EC numbers and descriptions as needed
    }

    formatted_hits = dbcan_hmmscan_formater(pd.read_csv(args.hits_csv), args.db_name, db_handler)
    formatted_hits.to_csv(args.output, index=False)
