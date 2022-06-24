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

CITATION = "None"
DOWNLOAD_OPTIONS ={'sulphur_tar_gz': {'version': VERSION}}
PROCESS_OPTIONS ={'sulphur_tar_gz': {'version': VERSION}}
DRAM_SETTINGS = {'sulphur_hmm': {'name': 'Sulphur placeholder', 'citation': CITATION,
                                 'notes': "This is just fegenie_hmm pretending to be sulphur."},
                  'sulphur_cutoffs': {'name': 'Sulphur cutoffs', 'citation': CITATION}
                 }

def download(temporary, logger, version=VERSION, verbose=True):
    """
    Retrieve genie release tar.gz

    This will get a tar file from the specified FeGenie release on git hub.

    :param temporary: Usually in the output dir
    :param verbose: TODO replace with logging setting
    :returns: Path to tar
    """
    logger.warn('This is not real and will just setup fegenie with a new name')
    NAME = 'FeGenie'
    database = path.join(temporary, f"{NAME}_{version}.tar.gz")
    # Note the 'v' in the name, GitHub wants it in the tag then it just takes it out. This could be a problem
    download_file(f"https://github.com/Arkadiy-Garber/FeGenie/archive/refs/tags/v{version}.tar.gz", logger,
                  database, verbose=verbose)
    return database


def process(input_file, output_dir, logger, threads=1,  version=VERSION, verbose=False) -> dict:
    temp_dir = path.dirname(input_file)
    # this is the path within the tar file
    NAME = 'FeGenie'
    tar_paths ={
        "sulphur_hmm":     [path.join(f"{NAME}-{version}", "iron", "iron_oxidation"), 
                            path.join(f"{NAME}-{version}", "iron", "iron_reduction")],
        "sulphur_cutoffs": path.join(f"{NAME}-{version}", "iron", "HMM-bitcutoffs.txt")
    }
    final_paths ={
        "sulphur_hmm"      : path.join(output_dir, f"{NAME}-{version}", "fegenie_iron_oxidation_reduction.hmm"),
        "sulphur_cutoffs" : path.join(output_dir, f"{NAME}-{version}", "fegenie_iron_cut_offs.txt")
    }

    new_fa_db_loc = path.join(output_dir, f"{NAME}_blast.faa")
    new_hmm_loc = path.join(output_dir, f"{NAME}_hmm.hmm")
    with tarfile.open(input_file, ) as tar:
        tar.extract(tar_paths["sulphur_cutoffs"], temp_dir)
        for info in tar.getmembers():
            tid = info.name
            if any([tid.startswith(i) for i in  tar_paths["sulphur_hmm"]]) and tid.endswith('hmm'):
                tar.extract(tid, temp_dir)
    
    # move and concatanate hmm to location
    if not path.exists(path.dirname(final_paths['sulphur_hmm'])):
        mkdir(path.dirname(final_paths['sulphur_hmm']))

    hmm_paths = [i for j in  tar_paths['sulphur_hmm'] for i in glob(path.join(temp_dir, j, '*.hmm'))]
    hmm_names = set() 
    with open(final_paths['sulphur_hmm'], 'wb') as wfd:
        for f in hmm_paths:
            if path.basename(f) not in hmm_names:
                hmm_names.add(path.basename(f))
                with open(f, 'rb') as fd:
                    copyfileobj(fd, wfd)

    # move the cutoffs
    move(path.join(temp_dir, tar_paths["sulphur_cutoffs"]), final_paths["sulphur_cutoffs"])
    
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

def hmmscan_formater(hits:pd.DataFrame,  db_name:str, hmm_info_path:str=None, top_hit:bool=True):
    if hmm_info_path is None:
        hmm_info = None
        hits_sig = hits[hits.apply(get_sig_row, axis=1)]
    else:
        hmm_info = pd.read_csv(hmm_info_path, sep='\t', index_col=0)
        hits_sig = sig_scores(hits, hmm_info)
    if len(hits_sig) == 0:
        # if nothing significant then return nothing, don't get descriptions
        return pd.DataFrame()
    if top_hit:
        # Get the best hits
        hits_sig = hits_sig.sort_values('full_evalue').drop_duplicates(subset=["query_id"])
    hits_df = hits_sig[['target_id', 'query_id', 'description']]
    hits_df.set_index('query_id', inplace=True, drop=True)
    hits_df.rename_axis(None, inplace=True)
    hits_df.columns = [f"{db_name}_id", f"{db_name}_description"]
    return hits_df


def search(genes_faa:str, tmp_dir:str, sulphur_hmm:str, sulphur_cutoffs:str, 
           logger:logging.Logger, threads:int, db_name:str=NAME, top_hit:bool=True, 
           verbose:bool=True):
    return run_hmmscan(genes_faa=genes_faa,
                       db_loc=sulphur_hmm,
                       db_name=db_name,
                       threads=threads,
                       output_loc=tmp_dir,
                       formater=partial(
                           hmmscan_formater,
                           db_name=db_name,
                           hmm_info_path=sulphur_cutoffs,
                           top_hit=True
                       ),
                       logger=logger)

