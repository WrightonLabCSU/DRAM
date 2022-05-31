from os import path
import tarfile
from shutil import move, rmtree
from mag_annotator.utils import download_file, run_process, make_mmseqs_db

CITATION = "CAMPER has no citeation and is in beta so you should not be using it."

DRAM_SETTINGS = { 
    'camper_hmm_loc':           {'origin': "camper_tar_gz", 'citation': CITATION, 'name': 'CAMPER HMM db'}, 
    'camper_fa_db_loc':         {'origin': "camper_tar_gz", 'citation': CITATION, 'name': 'CAMPER FASTA db'},
    'camper_hmm_cutoffs_loc':   {'origin': "camper_tar_gz", 'citation': CITATION, 'name': 'CAMPER HMM cutoffs'},
    'camper_distillate_loc':    {'origin': "camper_tar_gz", 'citation': CITATION, 'name': 'CAMPER Distillate form'},
    'camper_fa_db_cutoffs_loc': {'origin': "camper_tar_gz", 'citation': CITATION, 'name': 'CAMPER FASTA cutoffs'}}
DOWNLOAD_OPTIONS ={'version': '1.0.0-beta.1'}
PROCESS_OPTIONS ={'version': '1.0.0-beta.1'}
#NAME = 'CAMPER'

def download(temporary, logger, version=DOWNLOAD_OPTIONS['version'], verbose=True):
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


def process(camper_tar_gz_loc, output_dir, logger, version=PROCESS_OPTIONS['version'], 
                   threads=1, verbose=False) -> dict:
    name = f'CAMPER_{version}'
    temp_dir = path.dirname(camper_tar_gz_loc)
    tar_paths ={
        "camper_fa_db_loc"        : path.join(f"CAMPER-{version}", "CAMPER_blast.faa"),
        "camper_hmm_loc"          : path.join(f"CAMPER-{version}", "CAMPER.hmm"),
        "camper_fa_db_loc_scores" : path.join(f"CAMPER-{version}", "CAMPER_blast_scores.tsv"),
        "camper_distillate"       : path.join(f"CAMPER-{version}", "CAMPER_distillate.tsv"),
        "camper_hmm_cutoffs_loc"  : path.join(f"CAMPER-{version}", "CAMPER_hmm_scores.tsv"),
    }
    
    final_paths ={
        "camper_fa_db_loc"        : path.join(output_dir, "CAMPER_blast.faa"),
        "camper_hmm_loc"          : path.join(output_dir, "CAMPER.hmm"),
        "camper_fa_db_loc_scores" : path.join(output_dir, "CAMPER_blast_scores.tsv"),
        "camper_distillate"       : path.join(output_dir, "CAMPER_distillate.tsv"),
        "camper_hmm_cutoffs_loc"  : path.join(output_dir, "CAMPER_hmm_scores.tsv"),
    }
    
    new_fa_db_loc = path.join(output_dir, f"{name}_blast.faa")
    new_hmm_loc = path.join(output_dir, f"{name}_hmm.hmm")
    with tarfile.open(camper_tar_gz_loc) as tar:
        for v in tar_paths.values():
            tar.extract(v, temp_dir)
    
    # move tsv files, and hmm to location
    for i in ["camper_fa_db_loc_scores", "camper_distillate", "camper_hmm_cutoffs_loc", "camper_hmm_loc"]:
        move(path.join(temp_dir, tar_paths[i]), final_paths[i])
    
    # build dbs
    make_mmseqs_db(path.join(temp_dir, tar_paths["camper_fa_db_loc"]), final_paths["camper_fa_db_loc"], logger, threads=threads, verbose=verbose)
    run_process(['hmmpress', '-f', final_paths["camper_hmm_loc"]], logger, verbose=verbose)  # all are pressed just in case
    return {'fegenie_hmm': final_paths}

