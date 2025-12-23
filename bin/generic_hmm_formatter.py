#!/usr/bin/env python
import argparse
import pandas as pd
import click
import re

def extract_ec_numbers(definition):
    """
    Extract EC numbers from the definition string and format them into a semi-colon separated list.
    Each EC number is prefixed with 'EC:'.
    """
    # This regex pattern looks for EC numbers within square brackets and captures the numbers following "EC:"
    ec_numbers = re.findall(r'\[EC:(.*?)\]', definition)
    # Join the EC numbers with semi-colon separator and prepend 'EC:' to each number
    formatted_ec_numbers = '; '.join([f"EC:{ec.strip()}" for ec_block in ec_numbers for ec in ec_block.split()])
    return formatted_ec_numbers

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


@click.command()
@click.option("--hits_csv", type=click.Path(exists=True), help="Path to the HMM search results CSV file.")
@click.option("--hmm_info_path", type=click.Path(exists=True), help="Path to the hmm_info_path file.")
@click.option("--extract_ec_numbers", is_flag=True, help="Extract EC numbers from definitions. Can only be used if hmm_info_path is provided.")
@click.option("--gene_locs", type=click.Path(exists=True), help="Path to the gene locations TSV file.")
@click.option('--db_name', type=str, help='Name of the HMM database.')
@click.option("--output", type=click.Path(), help="Path to the formatted output file.")
def main(hits_csv, hmm_info_path, extract_ec_numbers, gene_locs, db_name, output):
# def generic_hmm_formatter(hits_csv, hmm_info_path, top_hit, output):
    hits_df = pd.read_csv(hits_csv)
    
    # Load hmm_info_path file and check if it's not empty (dummy sheet handling)
    if hmm_info_path is not None and not (hmm_info := pd.read_csv(hmm_info_path, sep="\t", index_col=0)).empty:
        # hmm_info = pd.read_csv(hmm_info_path, sep='\t', index_col=0)
        if extract_ec_numbers:
            hmm_info['definition'] = hmm_info['definition'].apply(extract_ec_numbers)
        hits_df = hits_df.merge(hmm_info, how='left', left_on="target_id", right_index=True)
        hits_df['bitScore'] = hits_df.apply(bitScore_per_row, axis=1)
        hits_df['score_rank'] = hits_df.apply(rank_per_row, axis=1)
        hits_df.dropna(subset=['score_rank'], inplace=True)
    # if not hmm_info_path/hmm_info_path is dummy sheet and extract_ec_numbers => error
    elif extract_ec_numbers:
        raise ValueError("extract_ec_numbers can only be used if hmm_info_path is provided.")
    else:

        hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)
        hits_df['score_rank'] = hits_df.apply(calculate_rank, axis=1)
        hits_df.dropna(subset=['score_rank'], inplace=True)

    hits_df['target_id'] = hits_df['target_id'].str.replace(r'.hmm', '', regex=True)
    
    if top_hit:
        hits_df = hits_df.sort_values('full_evalue').drop_duplicates(subset=["query_id"])
    
    hits_df.set_index('query_id', inplace=True, drop=True)
    hits_df.rename_axis(None, inplace=True)
    
    if 'definition' in hits_df.columns:
        hits_df = hits_df[['target_id', 'score_rank', 'bitScore', 'definition']]
        hits_df.columns = ["db_name_id", "db_name_rank", "db_name_bitScore", "db_name_hits"]
    else:
        hits_df = hits_df[['target_id', 'score_rank', 'bitScore']]
        hits_df.columns = ["db_name_id", "db_name_rank", "db_name_bitScore"]
    
    hits_df["db_name_search_type"] = 'hmm'
    
    hits_df.to_csv(output, index=True)

if __name__ == "__main__":
    main()