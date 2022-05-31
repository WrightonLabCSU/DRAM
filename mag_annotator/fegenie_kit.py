from os import path, mkdir
from glob import glob
import tarfile
from numpy import any
from shutil import rmtree, copyfileobj, move
from itertools import count
from mag_annotator.utils import download_file, run_process, make_mmseqs_db

CITATION = "Garber AI, Nealson KH, Okamoto A, McAllister SM, Chan CS, Barco RA and Merino N (2020) FeGenie: A Comprehensive Tool for the Identification of Iron Genes and Iron Gene Neighborhoods in Genome and Metagenome Assemblies. Front. Microbiol. 11:37. doi: 10.3389/fmicb.2020.00037"
NAME = 'FeGenie'
DOWNLOAD_OPTIONS ={'version': '1.2'}
PROCESS_OPTIONS ={'version': '1.2'}
DRAM_SETTINGS = {'fegenie_hmm': {'name': 'FeGenie', 'origin': "fegenie_tar_gz", 'citation': CITATION}}

def download(temporary, logger, version=DOWNLOAD_OPTIONS['version'], verbose=True):
    """
    Retrieve genie release tar.gz

    This will get a tar file from the specified FeGenie release on git hub.

    :param temporary: Usually in the output dir
    :param verbose: TODO replace with logging setting
    :returns: Path to tar
    """
    database = path.join(temporary, f"{NAME}_{version}.tar.gz")
    # Note the 'v' in the name, GitHub wants it in the tag then it just takes it out. This could be a problem
    download_file(f"https://github.com/Arkadiy-Garber/FeGenie/archive/refs/tags/v{version}.tar.gz", logger,
                  database, verbose=verbose)
    return database


def process(input_file, output_dir, logger, threads=1,  version=PROCESS_OPTIONS['version'], verbose=False) -> dict:
    temp_dir = path.dirname(input_file)
    # this is the path within the tar file
    tar_paths ={
        "fegenie_hmm":     [path.join(f"{NAME}-{version}", "iron", "iron_oxidation"), 
                            path.join(f"{NAME}-{version}", "iron", "iron_reduction")],
        "fegenie_cut_offs": path.join(f"{NAME}-{version}", "iron", "HMM-bitcutoffs.txt")
    }
    final_paths ={
        "fegenie_hmm"      : path.join(output_dir, f"{NAME}-{version}", "fegenie_iron_oxidation_reduction.hmm"),
        "fegenie_cut_offs" : path.join(output_dir, f"{NAME}-{version}", "fegenie_iron_cut_offs.txt")
    }

    new_fa_db_loc = path.join(output_dir, f"{NAME}_blast.faa")
    new_hmm_loc = path.join(output_dir, f"{NAME}_hmm.hmm")
    with tarfile.open(input_file, ) as tar:
        tar.extract(tar_paths["fegenie_cut_offs"], temp_dir)
        for info in tar.getmembers():
            tid = info.name
            if any([tid.startswith(i) for i in  tar_paths["fegenie_hmm"]]) and tid.endswith('hmm'):
                tar.extract(tid, temp_dir)
    
    # move and concatanate hmm to location
    if not path.exists(path.dirname(final_paths['fegenie_hmm'])):
        mkdir(path.dirname(final_paths['fegenie_hmm']))

    hmm_paths = [i for j in  tar_paths['fegenie_hmm'] for i in glob(path.join(temp_dir, j, '*.hmm'))]
    hmm_names = set() 
    with open(final_paths['fegenie_hmm'], 'wb') as wfd:
        for f in hmm_paths:
            if path.basename(f) not in hmm_names:
                hmm_names.add(path.basename(f))
                with open(f, 'rb') as fd:
                    copyfileobj(fd, wfd)
    breakpoint()

    # move the cutoffs
    move(path.join(temp_dir, tar_paths["fegenie_cut_offs"]), final_paths["fegenie_cut_offs"])
    
    # build dbs
    run_process(['hmmpress', '-f', final_paths["fegenie_hmm"]], logger, verbose=verbose)  # all are pressed just in case
    return final_paths

