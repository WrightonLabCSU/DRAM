from os import path, mkdir
from glob import glob
import tarfile
import logging
import pandas as pd
from numpy import any
from functools import partial
from shutil import rmtree, copyfileobj, move
from itertools import count
from mag_annotator.utils import download_file, run_process, make_mmseqs_db, \
    run_hmmscan, get_sig_row


VERSION = '1.2'
NAME = 'Sulphur'

CITATION = "Li W, O'Neill KR, Haft DH, DiCuccio M, Chetvernin V, Badretdin A, Coulouris G, Chitsaz F, Derbyshire MK, Durkin AS, Gonzales NR, Gwadz M, Lanczycki CJ, Song JS, Thanki N, Wang J, Yamashita RA, Yang M, Zheng C, Marchler-Bauer A, Thibaud-Nissen F. RefSeq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family model curation. Nucleic Acids Res. 2021 Jan 8;49(D1):D1020-D1028. doi: 10.1093/nar/gkaa1105. PMID: 33270901; PMCID: PMC7779008."
DOWNLOAD_OPTIONS ={'sulphur_tar_gz': {'version': VERSION}}
PROCESS_OPTIONS ={'sulphur_tar_gz': {'version': VERSION}}
DRAM_SETTINGS = {'sulphur_hmm': {'name': 'Sulphur placeholder', 'citation': CITATION,
                                 'notes': "This is just fegenie_hmm pretending to be sulphur."},
                 }

def process(input_file, output_dir, logger, threads=1,  version=VERSION, verbose=False) -> dict:
    # this is the path within the tar file
    final_paths ={
        "sulphur_hmm"      : path.join(output_dir, 
                                       f"{NAME}-{version}", f"{name}_hmm.hmm"),
    }
    if not path.exists(path.dirname(final_paths['sulphur_hmm'])):
        mkdir(path.dirname(final_paths['sulphur_hmm']))
    # move and concatanate hmm to location
    move(input_file, final_paths["sulphur_hmm"])
    
    # build dbs
    run_process(['hmmpress', '-f', final_paths["sulphur_hmm"]], logger, verbose=verbose)  # all are pressed just in case
    return final_paths

# TODO check this
def sig_scores(hits:pd.DataFrame, score_db:pd.DataFrame) -> pd.DataFrame:
    """
    This is a custom sig_scores function for FeGenie, it usese soft_bitscore_cutoff
    as a bit score cutoffs, given the name I am not shure that is corect.
    
    Also, I use full score, is that corect?
    """
    data = pd.merge(hits, score_db, how='left', left_on='target_id', right_index=True)
    return data[data['full_score'] > data['soft_bitscore_cutoff']]

def hmmscan_formater(hits:pd.DataFrame, logger:logging.Logger, db_name:str, top_hit:bool=True):
    """This is a formater for the results of a search"""
    hmm_info = None
    hits_sig = hits[hits.apply(get_sig_row, axis=1)]
    logger.debug(f"For {NAME}: there were {len(hits)} hits and {len(hits_sig)} were significant")
    if len(hits_sig) == 0:
        # if nothing significant then return nothing, don't get descriptions
        return pd.DataFrame(columns=[f"{db_name}_id"])
    if top_hit:
        # Get the best hits
        hits_sig = hits_sig.sort_values('full_evalue').drop_duplicates(subset=["query_id"])
    hits_df = hits_sig[['target_id', 'query_id']]
    hits_df.set_index('query_id', inplace=True, drop=True)
    hits_df.rename_axis(None, inplace=True)
    hits_df.columns = [f"{db_name}_id"]
    return hits_df


def search(genes_faa:str, tmp_dir:str, sulphur_hmm:str,  
           logger:logging.Logger, threads:int, db_name:str=NAME, top_hit:bool=True, 
           verbose:bool=True):
  return run_hmmscan(genes_faa=gene_faa,
                     db_loc=sulphur_hmm,
                     db_name=NAME,
                     threads=threads,
                     output_loc=tmp_dir,
                     formater=partial(
                         hmmscan_formater,
                         logger=logger,
                         db_name=NAME,
                         top_hit=True
                     ),
                     logger=logger)
    
# def search(query_db:str, gene_faa:str, tmp_dir:str, logger:logging.Logger, 
#            threads:str, verbose:str, db_handler, **args):
#     return run_hmmscan(genes_faa=gene_faa,
#                        db_loc=db_handler.config["search_databases"]['sulphur_hmm']['location'],
#                        db_name=NAME,
#                        threads=threads,
#                        output_loc=tmp_dir,
#                        formater=partial(
#                            hmmscan_formater,
#                            logger=logger,
#                            db_name=NAME,
#                            top_hit=True
#                        ),
#                        logger=logger)
