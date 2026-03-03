#!/usr/bin/env python
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from skbio.io import read as read_sequence
import os
from pathlib import Path
from utils.logger import get_logger
import click
from utils.click_utils import validate_comma_separated

FASTA_COLUMN = os.getenv('FASTA_COLUMN', 'input_fasta')

logger = get_logger(filename=Path(__file__).stem)

def read_and_preprocess(path: Path):
    input_fasta = input_fasta_from_filepath(path)
    try:
        df = pd.read_csv(path)
        df[FASTA_COLUMN] = input_fasta 
        return df
    except Exception as e:
        logger.error(f"Error loading DataFrame for input_fasta {input_fasta}: {str(e)}")
        return pd.DataFrame()

def input_fasta_from_filepath(file_path: Path):
    return file_path.stem.split("___")[0]

def assign_rank(row):
    rank = 'E'
    if row.get('kegg_bitScore', 0) > 350:
        rank = 'A'
    elif row.get('uniref_bitScore', 0) > 350:
        rank = 'B'
    elif row.get('kegg_bitScore', 0) > 60 or row.get('uniref_bitScore', 0) > 60:
        rank = 'C'
    elif any(row.get(f"{db}_bitScore", 0) > 60 for db in ['pfam', 'dbcan', 'merops']):
        rank = 'D'
    return rank

def convert_bit_scores_to_numeric(df):
    for col in df.columns:
        if "_bitScore" in col:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def count_motifs(gene_faa, motif="(C..CH)", genes_faa_dict=None):
    if genes_faa_dict is None:
        genes_faa_dict = dict()
    for seq in read_sequence(gene_faa, format="fasta"):
        if seq.metadata["id"] not in genes_faa_dict:
            genes_faa_dict[seq.metadata["id"]] = {}
        genes_faa_dict[seq.metadata["id"]]["heme_regulatory_motif_count"] = len(list(seq.find_with_regex(motif)))
    return genes_faa_dict

def set_gene_data(gene_faa, genes_faa_dict=None):
    if genes_faa_dict is None:
        genes_faa_dict = dict()
    for seq in read_sequence(gene_faa, format="fasta"):
        if seq.metadata["id"] not in genes_faa_dict:
            genes_faa_dict[seq.metadata["id"]] = {}
        split_label = seq.metadata["id"].split("_")
        gene_position = split_label[-1]
        start_position, end_position, strandedness = seq.metadata["description"].split("#")[1:4]

        input_fasta_name = Path(gene_faa).stem.split('_called_genes')[0]
        genes_faa_dict[seq.metadata["id"]][FASTA_COLUMN] = input_fasta_name
        genes_faa_dict[seq.metadata["id"]]["scaffold"] = (
            seq.metadata["id"]
            .removeprefix(genes_faa_dict[seq.metadata["id"]][FASTA_COLUMN])
            .removeprefix("_")
            .removesuffix(f"_{gene_position}"))
        genes_faa_dict[seq.metadata["id"]]["gene_number"] = int(gene_position)
        genes_faa_dict[seq.metadata["id"]]["start_position"] = int(start_position)
        genes_faa_dict[seq.metadata["id"]]["stop_position"] = int(end_position)
        genes_faa_dict[seq.metadata["id"]]["strandedness"] = int(strandedness)
    return genes_faa_dict

def organize_columns(df, special_columns=None):
    if special_columns is None:
        special_columns = []
    base_columns = ['query_id', FASTA_COLUMN, "scaffold",  'gene_number', 'start_position', 'stop_position', 'strandedness', 'rank']
    base_columns = [col for col in base_columns if col in df.columns]
    kegg_columns = sorted([col for col in df.columns if col.startswith('kegg_')], key=lambda x: (x != 'kegg_id', x))
    other_columns = [col for col in df.columns if col not in base_columns + kegg_columns + special_columns]
    db_prefixes = set(col.split('_')[0] for col in other_columns)
    sorted_other_columns = []
    for prefix in db_prefixes:
        prefixed_columns = sorted([col for col in other_columns if col.startswith(prefix + '_')], key=lambda x: (x != f"{prefix}_id", x))
        sorted_other_columns.extend(prefixed_columns)
    final_columns_order = base_columns + kegg_columns + sorted_other_columns + special_columns
    return df[final_columns_order]

@click.command()
@click.option("--annotations_dir", required=True, help="Directory of annotation files")
@click.option("--genes_dir", required=True, help="Directory genes faa file paths from prodigal")
@click.option("--output", help="Output file path for the combined annotations.")
@click.option("--threads", help="Number of threads for parallel processing", type=int, default=4)
def combine_annotations(annotations_dir, genes_dir, output, threads):
    """Combine annotation files with ranks and avoid duplicating specific columns."""
    annotations = Path(annotations_dir).glob("*")
    genes_faa = Path(genes_dir).glob("*")
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(read_and_preprocess, Path(path)) for path in annotations]
        data_frames = [future.result() for future in as_completed(futures)]
    
    combined_data = pd.concat([df for df in data_frames if not df.empty], ignore_index=True)
    if genes_faa:
        genes_faa_dict = dict()
        for gene_path in genes_faa:
            gene_path = str(gene_path)
            genes_faa_dict
            count_motifs(gene_path, "(C..CH)", genes_faa_dict=genes_faa_dict)
            set_gene_data(gene_path, genes_faa_dict)
        df = pd.DataFrame.from_dict(genes_faa_dict, orient='index').reset_index().rename(columns={'index': 'query_id'})
        combined_data = combined_data.drop(columns=df.columns.difference(["query_id", "scaffold", FASTA_COLUMN]), errors='ignore')
        combined_data = pd.merge(combined_data, df, how="outer", on=["query_id", FASTA_COLUMN])
    
    combined_data = convert_bit_scores_to_numeric(combined_data)

    aggregation_functions = {col: 'first' for col in combined_data.columns if col not in ['query_id', FASTA_COLUMN]}
    for col in ['Completeness', 'Contamination', 'taxonomy']:
        if col in combined_data.columns:
            aggregation_functions[col] = 'max'
    combined_data = combined_data.groupby(['query_id', FASTA_COLUMN], as_index=False).agg(aggregation_functions)
    # After aggregating data
    combined_data['rank'] = combined_data.apply(assign_rank, axis=1)

    # Continue with organizing columns and saving the DataFrame
    special_columns = ['Completeness', 'Contamination', 'taxonomy']
    special_columns = [col for col in special_columns if col in combined_data.columns]
    combined_data = organize_columns(combined_data, special_columns=special_columns)
    combined_data = combined_data.sort_values(by=[FASTA_COLUMN, 'scaffold', 'gene_number'])

    combined_data.to_csv(output, index=False, sep='\t')
    logger.info(f"Combined annotations saved to {output}, with corrected gene numbers.")

if __name__ == "__main__":
    combine_annotations()
