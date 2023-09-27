from os import path, stat
import tarfile
from shutil import move, rmtree
from mag_annotator.utils import download_file, run_process, make_mmseqs_db, \
    run_hmmscan, get_best_hits, BOUTFMT6_COLUMNS
from functools import partial
import logging
import pandas as pd

VERSION = '1.0.0-beta.1'
NAME = 'camper'
CITATION = "CAMPER has no citeation and is in beta so you should not be using it."
DRAM_SETTINGS = { 
    'camper_hmm':           {'citation': CITATION, 'name': 'CAMPER HMM db'}, 
    'camper_fa_db':         {'citation': CITATION, 'name': 'CAMPER FASTA db'},
    'camper_hmm_cutoffs':   {'citation': CITATION, 'name': 'CAMPER HMM cutoffs'},
    'camper_distillate':    {'citation': CITATION, 'name': 'CAMPER Distillate form'},
    'camper_fa_db_cutoffs': {'citation': CITATION, 'name': 'CAMPER FASTA cutoffs'}}
# the format is input file: options
DOWNLOAD_OPTIONS = {'camper_tar_gz': {'version': VERSION}}
PROCESS_OPTIONS = {'camper_tar_gz': {'version': VERSION}}
#NAME = 'CAMPER'

def download(temporary, logger, version=VERSION, verbose=True):
    """
    Retrieve CAMPER release tar.gz

    This will get a tar file that is automatically generated from making a campers release on git hub.  In order to 
    avoid changes in CAMPER being blindly excepted into DRAM, a new number must be put into the OPTIONS global
    variable in order to change this.

    :param temporary: Usually in the output dir
    :param verbose: TODO replace with logging setting
    :returns: Path to tar
    """
    camper_database = path.join(temporary, f"CAMPER_{version}.tar.gz")
    # Note the 'v' in the name, GitHub wants it in the tag then it just takes it out. This could be a problem
    download_file(f"https://github.com/WrightonLabCSU/CAMPER/archive/refs/tags/v{version}.tar.gz", logger,
                  camper_database, verbose=verbose)
    return camper_database


def process(camper_tar_gz, output_dir, logger, version=VERSION, 
                   threads=1, verbose=False) -> dict:
    name = f'CAMPER_{version}'
    temp_dir = path.dirname(camper_tar_gz)
    tar_paths ={
        "camper_fa_db"        : path.join(f"CAMPER-{version}", "CAMPER_blast.faa"),
        "camper_hmm"          : path.join(f"CAMPER-{version}", "CAMPER.hmm"),
        "camper_fa_db_cutoffs" : path.join(f"CAMPER-{version}", "CAMPER_blast_scores.tsv"),
        "camper_distillate"       : path.join(f"CAMPER-{version}", "CAMPER_distillate.tsv"),
        "camper_hmm_cutoffs"  : path.join(f"CAMPER-{version}", "CAMPER_hmm_scores.tsv"),
    }
    
    final_paths ={
        "camper_fa_db"        : path.join(output_dir, "CAMPER_blast.faa"),
        "camper_hmm"          : path.join(output_dir, "CAMPER.hmm"),
        "camper_fa_db_cutoffs" : path.join(output_dir, "CAMPER_blast_scores.tsv"),
        "camper_distillate"       : path.join(output_dir, "CAMPER_distillate.tsv"),
        "camper_hmm_cutoffs"  : path.join(output_dir, "CAMPER_hmm_scores.tsv"),
    }
    
    new_fa_db = path.join(output_dir, f"{name}_blast.faa")
    new_hmm = path.join(output_dir, f"{name}_hmm.hmm")
    with tarfile.open(camper_tar_gz) as tar:
        for v in tar_paths.values():
            tar.extract(v, temp_dir)
    
    # move tsv files, and hmm to location
    for i in ["camper_fa_db_cutoffs", "camper_distillate", "camper_hmm_cutoffs", "camper_hmm"]:
        move(path.join(temp_dir, tar_paths[i]), final_paths[i])
    
    # build dbs
    make_mmseqs_db(path.join(temp_dir, tar_paths["camper_fa_db"]), final_paths["camper_fa_db"], logger, threads=threads, verbose=verbose)
    run_process(['hmmpress', '-f', final_paths["camper_hmm"]], logger, verbose=verbose)  # all are pressed just in case
    return final_paths


def rank_per_row(row):
    r_a = row['A_rank']
    r_b = row['B_rank']
    score = row['bitScore']
    if score is None:
        return None
    if float(score) >= float(r_a):
         return 'A'
    if pd.isnull(r_b):
        return None
    if float(score) >= float(r_b):
         return 'B'
    return None


def blast_search_formater(hits_path, db_name, info_db, logger):
    if stat(hits_path).st_size == 0:
        return pd.DataFrame()
    hits = pd.read_csv(hits_path, sep='\t', header=None, 
                       names=BOUTFMT6_COLUMNS, index_col='qId')
    hits = hits.merge(info_db, how='left',left_on="tId", right_index=True)
    rank_col = f"{db_name}_rank"
    hits[rank_col] = hits.apply(rank_per_row, axis=1)
    hits.dropna(subset=[rank_col], inplace=True)
    logger.info('Getting descriptions of hits from %s' % db_name)
    hits = hits[['tId', rank_col, 'bitScore', 'ID_for_distillate', 'definition']]
    hits.rename(columns={
            'tId': f"{db_name}_hits",
            'ID_for_distillate': f"{db_name}_id", 
            'bitScore': f"{db_name}_bitScore",  
            'definition': f"{db_name}_definition"
        },
        inplace=True)
    hits[f"{db_name}_search_type"] = 'blast'
    return hits


def bitScore_per_row(row):
    if row['score_type'] == 'domain':
        return row.domain_score
    elif row['score_type'] == 'full':
        return row.full_score
    elif row['score_type'] == '-':
        return None
    else:
        raise ValueError("The score_type must be 'domain', 'full', or 'j")

#TODO decide if we need use_hmmer_thresholds:bool=False
def hmmscan_formater(hits:pd.DataFrame,  db_name:str, 
                             hmm_info_path:str, top_hit:bool=True):
    if hmm_info_path is None:
        hmm_info = None
        hits = hits[hits.apply(get_sig_row, axis=1)]
    else:
        hmm_info = pd.read_csv(hmm_info_path, sep='\t', index_col=0)
        hits = hits.merge(hmm_info, how='left',left_on="target_id", right_index=True)
        hits['bitScore'] = hits.apply(bitScore_per_row, axis=1)
        hits['score_rank'] = hits.apply(rank_per_row, axis=1)
        hits.dropna(subset=['score_rank'], inplace=True)
    if len(hits) == 0:
        # if nothing significant then return nothing, don't get descriptions
        return pd.DataFrame()
    if top_hit:
        # Get the best hits
        # TODO check we want top hit
        hits = hits.sort_values('full_evalue').drop_duplicates(subset=["query_id"])
    hits.set_index('query_id', inplace=True, drop=True)
    hits.rename_axis(None, inplace=True)
    if 'definition' in hits.columns:
        hits = hits[['target_id', 'score_rank', 'bitScore', 'definition']]
        hits.columns = [f"{db_name}_id", f"{db_name}_rank", 
                        f"{db_name}_bitScore", f"{db_name}_hits"]
    else:
        hits = hits[['target_id', 'score_rank', 'bitScore']]
        hits.columns = [f"{db_name}_id", f"{db_name}_rank", 
                        f"{db_name}_bitScore"]
    # Rename
    hits[f"{db_name}_search_type"] = 'hmm'
    return hits


def get_minimum_bitscore(info_db):
    bit_score_threshold = min(info_db[['A_rank', 'B_rank']].min().values)
    return bit_score_threshold


def blast_search(query_db, target_db, working_dir, info_db_path, 
                 db_name, logger, threads=10, verbose=False):
    """A convenience function to do a blast style forward best hits search"""
    # Get kegg hits
    info_db = pd.read_csv(info_db_path, sep='\t', index_col=0)
    bit_score_threshold = get_minimum_bitscore(info_db)
    hits_path = get_best_hits(query_db, target_db, logger, working_dir, 'gene', db_name, 
                         bit_score_threshold, threads, verbose=verbose)
    return blast_search_formater(hits_path, db_name, info_db, logger)


# in the future the database will get the same input as was given in the data
def search(query_db:str, genes_faa:str, tmp_dir:str, logger:logging.Logger, 
           threads:str, verbose:str, camper_fa_db:str, camper_hmm:str, 
           camper_fa_db_cutoffs:str, camper_hmm_cutoffs:str):
        fasta = blast_search(query_db=query_db, 
                             target_db=camper_fa_db, 
                             working_dir=tmp_dir, 
                             info_db_path=camper_fa_db_cutoffs,
                             db_name=NAME, 
                             logger=logger,
                             threads=threads,
                             verbose=verbose)
        hmm = run_hmmscan(genes_faa=genes_faa,
                          db_loc=camper_hmm,
                          db_name=NAME,
                          threads=threads,
                          output_loc=tmp_dir,
                          logger=logger,
                          formater=partial(
                              hmmscan_formater,
                              db_name=NAME,
                              hmm_info_path=camper_hmm_cutoffs,
                              top_hit=True
                          ))
        full = pd.concat([fasta, hmm])
        if len(full) < 1:
            return pd.DataFrame()
        return (full
               .groupby(full.index)
               .apply(
                   lambda x: (x
                              .sort_values(f'{NAME}_search_type', ascending=True) # make sure hmm is first
                              .sort_values(f'{NAME}_bitScore', ascending=False)
                              .iloc[0])
                  )
               )

