"""
Removing duplicates and adding our own headers to files. Dram dose this with a genome link but I wanted an independent tool.

RE-HEADER

I don't know if kegg used to do the headers differently but that is what I 
assume is the case, because simply concatenating the pep files dose not work
anymore. 

For now it seems simplest to just use the information in the kff files to add 
the information back in the same format. I would like a beter long term solution
but this is the world we live in, so be it. 



FROM: https://www.kegg.jp/kegg/download/Readme/README.genes

The .kff file contains the following information for each gene in the genome:

   column 1  - KEGG GENES identifier (in the form of org:gene)
   column 2  - Feature key such as CDS, tRNA, rRNA, ncRNA, lncRNA, miRNA, and
               gene (usually meaning pseudogene)
   column 3  - Amino acid sequence length
   column 4  - Nucleotide sequence length
   column 5  - Geonomic position
   column 6  - NCBI GeneID
   column 7  - NCBI ProteinID
   column 8  - Gene name
   column 9  - Definition in the original DB
   column 10 - KO identifier (K number)
   column 11 - KO definition

In order for kegg to work it needs to be in the fallowing format

  column 1  - KEGG GENES identifier (in the form of org:gene)
    Falowed by a space
  column 2  - A kegg id in the format (k\d\d\d\d\d) parentheses included
  column 3  - KO definition with EC numbers that have format \[EC:\d*.\d*.\d*.\d*\]'
  then a ';' and any other information you want I don't think it maters. 
"""
import logging
from glob import glob
from multiprocessing import Pool
from os.path import join
from Bio import SeqIO
from functools import partial
import pandas as pd
import numpy as np
import argparse
# from ftplib import FTP
# from getpass import getpass

  

# wget -r ftp://$UNAME:$PASS@ftp.kegg.net/kegg/README.kegg
# wget -r ftp://$UNAME:$PASS@ftp.kegg.net/kegg/RELEASE
# wget -r ftp://$UNAME:$PASS@ftp.kegg.net/kegg/genes/MD5.genes
# wget -r  -A "*\.pep\.gz" ftp://$UNAME:$PASS@ftp.kegg.net/kegg/genes/organisms
# wget -r  -A "*\.kff\.gz" ftp://$UNAME:$PASS@ftp.kegg.net/kegg/genes/organisms
# 
# uname = input("What is the kegg User Name: ")
# passw = getpass("What is the kegg password: ")
# 
# with FTP('ftp.kegg.net', user=uname, passwd=passw) as ftp:
#     ftp.login()
#     # changing directory
#     ftp.retrlines('LIST')
# 
# # Done

class CallCounted:
    """Decorator to determine number of calls for a method

    nice steel https://stackoverflow.com/questions/812477/how-many-times-was-logging-error-called
    """

    def __init__(self,method):
        self.method=method
        self.counter=0

    def __call__(self,*args,**kwargs):
        self.counter+=1
        return self.method(*args,**kwargs)


logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
logging.error = CallCounted(logging.error)
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# set a format which is simpler for console use
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
# tell the handler to use this format
console.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger().addHandler(console)



def get_file_path(par_path:str, ext:str):
    if len(chi_paths := glob(join(par_path, f"*.{ext}"))) < 1:
        return None
    elif len(chi_paths) > 1:
        logging.error( "%s", f" The folder {par_path} multi {ext}")
        return None
    else:
        chi_path = chi_paths[0]
        return chi_path


def get_file_sets(organisms):
    for path in glob(join(organisms, "*")):
        kff, pep = (get_file_path(path, 'kff'), get_file_path(path, 'pep'))
        if kff is None and pep is None:
            logging.error("%s", f" The folder {path} has no data")
            continue
        if pep is None:
            logging.error("%s", f" The folder {path} has a kff but no pep files")
            continue
        if kff is None:
            logging.error("%s", f" The folder {path} has a pep but no kff files")
            continue
        yield (kff, pep)


def re_header(x, y):
    x.description = y
    return x


def check_kff_pep(kff:set, pep:set):
    if (uniques := len(kff - pep)) > 0:
        pass
        # logging.error("%s", f"{uniques} kegg ids are in the kff, but not in the pep")
    if (uniques := len(pep - kff)) > 0:
        logging.error("%s", f"{uniques} kegg ids are in the pep, but not in the kff")


def make_kff_args(kff_str:str):
    kff_line = [None, None, None]
    logging.info("%s", f"Processing kff line : {kff_str}")
    try:
        kff_line[0] = kff_str.split('\t')[0]
        kff_line[1] = kff_str.split('\t')[9]
        kff_line[2] = kff_str.split('\t')[10]
    except IndexError:
        logging.warn("%s", f"Possible illegitimate kff line : {kff_line}")
    return kff_line


def read_join_pep_kff(paths:tuple):
    kff_path, pep_path  = paths
    kff_rws = np.genfromtxt(kff_path, dtype='str', delimiter='\n')
    kff_arr = [make_kff_args(kff_str) for kff_str in kff_rws]
    kff_df = pd.DataFrame(
        kff_arr,
        columns=['kegg_id', 'ko_id', 'ko_description']
    )
    pep_seq = np.array([record for record in SeqIO.parse(pep_path, 'fasta')], dtype=object)
    pep_ids = np.array([record.id for record in pep_seq], dtype=object)
    pep_idf = pd.DataFrame(pep_ids, columns=['kegg_id'])
    check_kff_pep(set(kff_df['kegg_id']), set(pep_ids))
    headers = pd.merge(pep_idf, kff_df, how='left', on='kegg_id')
    headers = headers.apply(lambda x: f"{x['kegg_id']} ({x['ko_id']}) {x['ko_description']}", axis=1).values
    pep_out = [re_header(i, j) for i, j in zip(pep_seq, headers)]

    return np.array([pep_ids, pep_out], dtype=object)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Remove duplicates from KEGG pep files")
    parser.add_argument("--pep_loc", type=str, help="Path to the KEGG pep files parent dir", default=join("ftp.kegg.net", "kegg", "genes", "organisms"))
    parser.add_argument("--output_file", type=str, help="Path to the output file", default="kegg-all-orgs_unique_reheader.pep")
    args = parser.parse_args()

    organisms = [i for i in get_file_sets(args.pep_loc)]
    if logging.error.counter > 0:
        logging.critical("The program will end execution because there may be errors in downloading.")
        raise ValueError("System exiting see logs for why.")
    with Pool(30) as pool:
        pep = np.concatenate(pool.map(read_join_pep_kff, organisms), axis=1, dtype=object)

    ids=pep[0,:]
    unique_ids, indices = np.unique(ids, return_index=True)
    pep_out=pep[1,indices]
    logging.info("%s", f"The number of non unique organisms is: {len(ids) - len(pep_out)}")
    with open(args.output_file, 'a') as outFile:
        SeqIO.write(pep_out, outFile, 'fasta')
