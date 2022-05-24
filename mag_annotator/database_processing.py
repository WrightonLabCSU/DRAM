"""
Contains most of the backend for the DRAM_setup.py script, used to setup databases for each user.
"""
from os import path, mkdir
from datetime import datetime
from shutil import move, rmtree
from glob import glob
import gzip
import tarfile
from collections import defaultdict
import logging

from skbio import read as read_sequence
from skbio import write as write_sequence

from mag_annotator.database_handler import DatabaseHandler
from mag_annotator.utils import run_process, make_mmseqs_db, download_file, merge_files, remove_prefix, setup_logger


DEFAULT_DBCAN_RELEASE = '10'
DEFAULT_DBCAN_DATE = '07292021'
DEFAULT_UNIREF_VERSION = '90'
DFLT_OUTPUT_DIR = '.'
CAMPER_RELEASE = '1.0.0-beta.1'
LOGGER = logging.getLogger("database_processing.log")

# TODO: check if dbcan or pfam is down, raise appropriate error
# TODO: upgrade to pigz?


def get_iso_date():
    return datetime.today().strftime('%Y%m%d')


def download_pfam_descriptions(output_dir='.', verbose=True):
    pfam_hmm_dat = path.join(output_dir, 'Pfam-A.hmm.dat.gz')
    download_file('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz', pfam_hmm_dat,
                  verbose=verbose)
    return pfam_hmm_dat


def download_dbcan(dbcan_hmm=None, output_dir='.', dbcan_release=DEFAULT_DBCAN_RELEASE, verbose=True):
    dbcan_hmm = path.join(output_dir, f"dbCAN-HMMdb-V{dbcan_release}.txt" )
    if int(dbcan_release) < int(DEFAULT_DBCAN_RELEASE):
        link_path = f"http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V{dbcan_release}.txt"
    else:
        link_path = f"http://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V{dbcan_release}.txt"

    print(f"Downloading dbCAN from: {link_path}")
    download_file(link_path, dbcan_hmm, verbose=verbose)
    return dbcan_hmm


def download_dbcan_descriptions(output_dir='.', dbcan_release=DEFAULT_DBCAN_RELEASE, upload_date=DEFAULT_DBCAN_DATE, 
                                verbose=True):
    dbcan_fam_activities = path.join(output_dir, f'CAZyDB.{upload_date}.fam-activities.txt')
    link_path = f"https://bcb.unl.edu/dbCAN2/download/Databases/V{dbcan_release}/CAZyDB.{upload_date}.fam-activities.txt"
    print(f"Downloading dbCAN family activities from : {link_path}")
    download_file(link_path, dbcan_fam_activities, verbose=verbose)
    return dbcan_fam_activities


def download_dbcan_subfam_ec(output_dir='.', dbcan_release=DEFAULT_DBCAN_RELEASE, upload_date=DEFAULT_DBCAN_DATE, verbose=True):
    dbcan_subfam_ec = path.join(output_dir, f"CAZyDB.{upload_date}.fam.subfam.ec.txt")
    link_path = (f"https://bcb.unl.edu/dbCAN2/download/Databases/"
                 f"V{dbcan_release}/CAZyDB.{upload_date}.fam.subfam.ec.txt")
    print(f"Downloading dbCAN sub-family encumber from : {link_path}")
    download_file(link_path, dbcan_subfam_ec, verbose=verbose)
    return dbcan_subfam_ec


def download_kofam_hmms(output_dir='.', verbose=False):
    kofam_profile_tar_gz = path.join(output_dir, 'kofam_profiles.tar.gz')
    download_file('ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz', kofam_profile_tar_gz, verbose=verbose)
    return kofam_profile_tar_gz

def download_kofam_hmms(output_dir='.', verbose=False):
    kofam_profile_tar_gz = path.join(output_dir, 'kofam_profiles.tar.gz')
    download_file('ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz', kofam_profile_tar_gz, verbose=verbose)
    return kofam_profile_tar_gz

def generate_modified_kegg_fasta(kegg_fasta, gene_ko_link_loc=None):
    """Takes kegg fasta file and gene ko link file, adds kos not already in headers to headers"""
    genes_ko_dict = defaultdict(list)
    if gene_ko_link_loc is not None:
        if gene_ko_link_loc.endswith('.gz'):
            gene_ko_link_fh = gzip.open(gene_ko_link_loc, 'rt')
        else:
            gene_ko_link_fh = open(gene_ko_link_loc)
        for line in gene_ko_link_fh:
            gene, ko = line.strip().split()
            genes_ko_dict[gene].append(remove_prefix(ko, 'ko:'))
    for seq in read_sequence(kegg_fasta, format='fasta'):
        new_description = seq.metadata['description']
        for ko in genes_ko_dict[seq.metadata['id']]:
            if ko not in new_description:
                new_description += '; %s' % ko
        seq.metadata['description'] = new_description
        yield seq


def process_kegg_db(output_dir, kegg_loc, gene_ko_link_loc=None, download_date=None, threads=10, verbose=True):
    if download_date is None:
        download_date = get_iso_date()
    if gene_ko_link_loc is not None:
        # add KOs to end of header where KO is not already there
        kegg_mod_loc = path.join(output_dir, 'kegg.mod.fa')
        write_sequence(generate_modified_kegg_fasta(kegg_loc, gene_ko_link_loc),
                       format='fasta', into=kegg_mod_loc)
    else:
        kegg_mod_loc = kegg_loc
    # make mmseqsdb from modified kegg fasta
    kegg_mmseqs_db = path.join(output_dir, 'kegg.%s.mmsdb' % download_date)
    make_mmseqs_db(kegg_mod_loc, kegg_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
    return kegg_mmseqs_db

def process_camper(camper_tar_gz_loc, output_dir, 
                   threads=1, verbose=False) -> dict:
    name = f'CAMPER_{CAMPER_RELEASE}'
    temp_dir = path.dirname(camper_tar_gz_loc)
    tar_paths ={
        "camper_blast"       : path.join(f"CAMPER-{CAMPER_RELEASE}", "CAMPER_blast.faa"),
        "camper_hmm"             : path.join(f"CAMPER-{CAMPER_RELEASE}", "CAMPER.hmm"),
        "camper_blast_scores": path.join(f"CAMPER-{CAMPER_RELEASE}", "CAMPER_blast_scores.tsv"),
        "camper_distillate"  : path.join(f"CAMPER-{CAMPER_RELEASE}", "CAMPER_distillate.tsv"),
        "camper_hmm_scores"  : path.join(f"CAMPER-{CAMPER_RELEASE}", "CAMPER_hmm_scores.tsv"),
    }
    
    final_paths ={
        "camper_blast"       : path.join(output_dir, "CAMPER_blast.faa"),
        "camper_hmm"         : path.join(output_dir, "CAMPER.hmm"),
        "camper_blast_scores": path.join(output_dir, "CAMPER_blast_scores.tsv"),
        "camper_distillate"  : path.join(output_dir, "CAMPER_distillate.tsv"),
        "camper_hmm_scores"  : path.join(output_dir, "CAMPER_hmm_scores.tsv"),
    }
    
    new_fa_db_loc = path.join(output_dir, f"{name}_blast.faa")
    new_hmm_loc = path.join(output_dir, f"{name}_hmm.hmm")
    with tarfile.open(camper_tar_gz_loc) as tar:
        for v in tar_paths.values():
            tar.extract(v, temp_dir)
    
    # move tsv files, and hmm to location
    for i in ["camper_blast_scores", "camper_distillate", "camper_hmm_scores", "camper_hmm"]:
        move(path.join(temp_dir, tar_paths[i]), final_paths[i])
    
    # build dbs
    make_mmseqs_db(path.join(temp_dir, tar_paths["camper_blast"]), final_paths["camper_blast"], threads=threads, verbose=verbose)
    run_process(['hmmpress', '-f', final_paths["camper_hmm"]], verbose=verbose)  # all are pressed just in case
    return final_paths

def process_kofam_hmms(kofam_profile_tar_gz, output_dir=DFLT_OUTPUT_DIR, verbose=False):
    kofam_profiles = path.join(output_dir, 'kofam_profiles')
    mkdir(kofam_profiles)
    run_process(['tar', '-xzf', kofam_profile_tar_gz, '-C', kofam_profiles], verbose=verbose)
    merged_kofam_profiles = path.join(output_dir, 'kofam_profiles.hmm')
    merge_files(glob(path.join(kofam_profiles, 'profiles', '*.hmm')), merged_kofam_profiles)
    run_process(['hmmpress', '-f', merged_kofam_profiles], verbose=verbose)
    return merged_kofam_profiles


def download_kofam_ko_list(output_dir='.', verbose=False):
    kofam_ko_list_gz = path.join(output_dir, 'kofam_ko_list.tsv.gz')
    download_file('ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz', kofam_ko_list_gz, verbose=verbose)


def process_kofam_ko_list(kofam_ko_list_gz, output_dir='.', verbose=False):
    # TODO: fix this so that it is gunzipped to the path
    kofam_ko_list = path.join(output_dir, 'kofam_ko_list.tsv')
    run_process(['gunzip', '-c', kofam_ko_list_gz, '>', kofam_ko_list], verbose=verbose)
    return kofam_ko_list


# def download_and_process_kofam_ko_list(kofam_ko_list_gz=None, output_dir='.', verbose=False):
#     if kofam_ko_list_gz is None:
#         kofam_ko_list_gz = path.join(output_dir, 'kofam_ko_list.tsv.gz')
#         download_file('ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz', kofam_ko_list_gz, verbose=verbose)
#     # TODO: fix this so that it is gunzipped to the path
#     kofam_ko_list = path.join(output_dir, 'kofam_ko_list.tsv')
#     run_process(['gunzip', kofam_ko_list_gz], verbose=verbose)
#     return kofam_ko_list

def download_uniref(output_dir='.', uniref_version=DEFAULT_UNIREF_VERSION, 
                    threads=10, verbose=True):
    uniref_fasta_zipped = path.join(output_dir, 'uniref%s.fasta.gz' % uniref_version)
    uniref_url = 'https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz' % \
                 (uniref_version, uniref_version)
    download_file(uniref_url, uniref_fasta_zipped, verbose=verbose)
    return uniref_fasta_zipped

def process_uniref(uniref_fasta_zipped, output_dir='.', 
                   uniref_version=DEFAULT_UNIREF_VERSION, threads=10,
                   verbose=True):
    uniref_mmseqs_db = path.join(output_dir, 'uniref%s.%s.mmsdb' % (uniref_version, get_iso_date()))
    make_mmseqs_db(uniref_fasta_zipped, uniref_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
    return uniref_mmseqs_db

# def download_and_process_uniref(uniref_fasta_zipped=None, output_dir='.', 
#                                 uniref_version=DEFAULT_UNIREF_VERSION, threads=10,
#                                 verbose=True):
#     if uniref_fasta_zipped is None:  # download database if not provided
#         uniref_fasta_zipped = path.join(output_dir, 'uniref%s.fasta.gz' % uniref_version)
#         uniref_url = 'https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz' % \
#                      (uniref_version, uniref_version)
#         download_file(uniref_url, uniref_fasta_zipped, verbose=verbose)
#     uniref_mmseqs_db = path.join(output_dir, 'uniref%s.%s.mmsdb' % (uniref_version, get_iso_date()))
#     make_mmseqs_db(uniref_fasta_zipped, uniref_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
#     return uniref_mmseqs_db


def process_mmspro(full_alignment, output_dir, db_name='db', threads=10, verbose=True):
    mmseqs_msa = path.join(output_dir, '%s.mmsmsa' % db_name)
    run_process(['mmseqs', 'convertmsa', full_alignment, mmseqs_msa], verbose=verbose)
    mmseqs_profile = path.join(output_dir, '%s.mmspro' % db_name)
    run_process(['mmseqs', 'msa2profile', mmseqs_msa, mmseqs_profile, '--match-mode', '1', '--threads', str(threads)],
                verbose=verbose)
    tmp_dir = path.join(output_dir, 'tmp')
    run_process(['mmseqs', 'createindex', mmseqs_profile, tmp_dir, '-k', '5', '-s', '7', '--threads', str(threads)],
                verbose=verbose)
    return mmseqs_profile


def download_pfam(output_dir='.', verbose=True):
    pfam_full_zipped = path.join(output_dir, 'Pfam-A.full.gz')
    download_file('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz', pfam_full_zipped,
                  verbose=verbose)
    return pfam_full_zipped


def process_pfam(pfam_full_zipped, output_dir='.', threads=10, verbose=True):
    pfam_profile = process_mmspro(pfam_full_zipped, output_dir, 'pfam', threads, verbose)
    return pfam_profile


# def download_and_process_pfam(pfam_full_zipped=None, output_dir='.', threads=10, verbose=True):
#     if pfam_full_zipped is None:  # download database if not provided
#         pfam_full_zipped = path.join(output_dir, 'Pfam-A.full.gz')
#         download_file('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz', pfam_full_zipped,
#                       verbose=verbose)
#     pfam_profile = process_mmspro(pfam_full_zipped, output_dir, 'pfam', threads, verbose)
#     return pfam_profile


def process_dbcan(dbcan_hmm, verbose=True):
    run_process(['hmmpress', '-f', dbcan_hmm], verbose=verbose)
    return dbcan_hmm


def download_viral_refseq(output_dir='.', viral_files=2, verbose=True):
    """Can only download newest version"""
    # download all of the viral protein files, need to know the number of files
    # TODO: Make it so that you don't need to know number of viral files in refseq viral

    faa_base_name = 'viral.%s.protein.faa.gz'
    viral_faa_glob = path.join(output_dir, faa_base_name % '*')
    for number in range(viral_files):
        number += 1
        refseq_url = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.%s.protein.faa.gz' % number
        refseq_faa = path.join(output_dir, faa_base_name % number)
        download_file(refseq_url, refseq_faa, verbose=verbose)

    # then merge files from above
    merged_viral_faas = path.join(output_dir, 'viral.merged.protein.faa.gz')
    run_process(['cat %s > %s' % (' '.join(glob(viral_faa_glob)), merged_viral_faas)], shell=True)
    return merged_viral_faas


def process_viral_refseq(merged_viral_faas, output_dir='.', viral_files=2, threads=10, verbose=True):
    refseq_viral_mmseqs_db = path.join(output_dir, 'refseq_viral.%s.mmsdb' % get_iso_date())
    make_mmseqs_db(merged_viral_faas, refseq_viral_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
    return refseq_viral_mmseqs_db

# def download_and_process_viral_refseq(merged_viral_faas=None, output_dir='.', viral_files=2, threads=10, verbose=True):
#     """Can only download newest version"""
#     # download all of the viral protein files, need to know the number of files
#     # TODO: Make it so that you don't need to know number of viral files in refseq viral
# 
#     if merged_viral_faas is None:  # download database if not provided
#         faa_base_name = 'viral.%s.protein.faa.gz'
#         viral_faa_glob = path.join(output_dir, faa_base_name % '*')
#         for number in range(viral_files):
#             number += 1
#             refseq_url = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.%s.protein.faa.gz' % number
#             refseq_faa = path.join(output_dir, faa_base_name % number)
#             download_file(refseq_url, refseq_faa, verbose=verbose)
# 
#         # then merge files from above
#         merged_viral_faas = path.join(output_dir, 'viral.merged.protein.faa.gz')
#         run_process(['cat %s > %s' % (' '.join(glob(viral_faa_glob)), merged_viral_faas)], shell=True)
# 
#     # make mmseqs database
#     refseq_viral_mmseqs_db = path.join(output_dir, 'refseq_viral.%s.mmsdb' % get_iso_date())
#     make_mmseqs_db(merged_viral_faas, refseq_viral_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
#     return refseq_viral_mmseqs_db


def download_merops_peptidases(output_dir='.', verbose=True):
    peptidase_faa = path.join(output_dir, 'merops_peptidases_nr.faa')
    merops_url = 'ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib'
    download_file(merops_url, peptidase_faa, verbose=verbose)
    return peptidase_faa


def process_merops_peptidases(peptidase_faa, output_dir='.', threads=10, verbose=True):
    peptidase_mmseqs_db = path.join(output_dir, 'peptidases.%s.mmsdb' % get_iso_date())
    make_mmseqs_db(peptidase_faa, peptidase_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
    return peptidase_mmseqs_db


# def download_and_process_merops_peptidases(peptidase_faa=None, output_dir='.', threads=10, verbose=True):
#     if peptidase_faa is None:  # download database if not provided
#         peptidase_faa = path.join(output_dir, 'merops_peptidases_nr.faa')
#         merops_url = 'ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib'
#         download_file(merops_url, peptidase_faa, verbose=verbose)
#     peptidase_mmseqs_db = path.join(output_dir, 'peptidases.%s.mmsdb' % get_iso_date())
#     make_mmseqs_db(peptidase_faa, peptidase_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
#     return peptidase_mmseqs_db


def download_vogdb(output_dir='.', vogdb_release='latest', verbose=True):
    vog_hmm_targz = path.join(output_dir, 'vog.hmm.tar.gz')
    vogdb_url = 'http://fileshare.csb.univie.ac.at/vog/%s/vog.hmm.tar.gz' % vogdb_release
    download_file(vogdb_url, vog_hmm_targz, verbose=verbose)
    return vog_hmm_targz


def process_vogdb(vog_hmm_targz, output_dir='.', vogdb_release='latest', verbose=True):
    hmm_dir = path.join(output_dir, 'vogdb_hmms')
    mkdir(hmm_dir)
    vogdb_targz = tarfile.open(vog_hmm_targz)
    vogdb_targz.extractall(hmm_dir)
    vog_hmms = path.join(output_dir, 'vog_%s_hmms.txt' % vogdb_release)
    merge_files(glob(path.join(hmm_dir, 'VOG*.hmm')), vog_hmms)
    run_process(['hmmpress', '-f', vog_hmms], verbose=verbose)
    return vog_hmms


# def download_and_process_vogdb(vog_hmm_targz=None, output_dir='.', vogdb_release='latest', verbose=True):
#     if vog_hmm_targz is None:
#         vog_hmm_targz = path.join(output_dir, 'vog.hmm.tar.gz')
#         vogdb_url = 'http://fileshare.csb.univie.ac.at/vog/%s/vog.hmm.tar.gz' % vogdb_release
#         download_file(vogdb_url, vog_hmm_targz, verbose=verbose)
#     hmm_dir = path.join(output_dir, 'vogdb_hmms')
#     mkdir(hmm_dir)
#     vogdb_targz = tarfile.open(vog_hmm_targz)
#     vogdb_targz.extractall(hmm_dir)
#     vog_hmms = path.join(output_dir, 'vog_%s_hmms.txt' % vogdb_release)
#     merge_files(glob(path.join(hmm_dir, 'VOG*.hmm')), vog_hmms)
#     run_process(['hmmpress', '-f', vog_hmms], verbose=verbose)
#     return vog_hmms


def download_vog_annotations(output_dir, vogdb_version='latest', verbose=True):
    vog_annotations = path.join(output_dir, 'vog_annotations_%s.tsv.gz' % vogdb_version)
    download_file('http://fileshare.csb.univie.ac.at/vog/%s/vog.annotations.tsv.gz' % vogdb_version,
                  vog_annotations, verbose=verbose)
    return vog_annotations


def download_genome_summary_form(output_dir, branch='master', verbose=True):
    genome_summary_form = path.join(output_dir, 'genome_summary_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/genome_summary_form.tsv' % branch,
                  genome_summary_form, verbose=verbose)
    return genome_summary_form


def download_module_step_form(output_dir, branch='master', verbose=True):
    function_heatmap_form = path.join(output_dir, 'module_step_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/module_step_form.tsv' % branch,
                  function_heatmap_form, verbose=verbose)
    return function_heatmap_form


def download_etc_module_database(output_dir, branch='master', verbose=True):
    etc_module_database = path.join(output_dir, 'etc_mdoule_database.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/etc_module_database.tsv' % branch,
                  etc_module_database, verbose=verbose)
    return etc_module_database


def download_function_heatmap_form(output_dir, branch='master', verbose=True):
    function_heatmap_form = path.join(output_dir, 'function_heatmap_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/function_heatmap_form.tsv' % branch,
                  function_heatmap_form, verbose=verbose)
    return function_heatmap_form


def download_amg_database(output_dir, branch='master', verbose=True):
    amg_database = path.join(output_dir, 'amg_database.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/amg_database.tsv' % branch,
                  amg_database, verbose=verbose)
    return amg_database


def check_file_exists(*paths):
    for i in paths:
        if i is None:
            continue
        elif path.isfile(i):
            continue
        else:
            raise ValueError(f"Database location does not exist: {i}")


def download_camper_tar_gz(temporary, verbose=True):
    """
    Retrieve CAMPER release tar.gz

    This will get a tar file that is automatically generated from making a campers release on git hub.  In order to 
    avoid changes in CAMPER being blindly excepted into DRAM, a new number must be put into the CAMPER_RELEASE global
    variable in order to change this.

    :param temporary: Usually in the output dir
    :param verbose: TODO replace with logging setting
    :returns: Path to tar
    """
    camper_database = path.join(temporary, f"CAMPER_{CAMPER_RELEASE}.tar.gz")
    # Note the 'v' in the name, GitHub wants it in the tag then it just takes it out. This could be a problem
    download_file(f"https://github.com/WrightonLabCSU/CAMPER/archive/refs/tags/v{CAMPER_RELEASE}.tar.gz",
                  camper_database, verbose=verbose)
    return camper_database


def prepare_databases(output_dir, kegg_loc=None, gene_ko_link_loc=None, kofam_hmm_loc=None,
                      kofam_ko_list_loc=None, kegg_download_date=None, uniref_loc=None,
                      uniref_version=DEFAULT_UNIREF_VERSION, pfam_loc=None, pfam_hmm_dat=None,
                      dbcan_loc=None, dbcan_fam_activities:str=None, dbcan_subfam_ec:str=None, dbcan_version=DEFAULT_DBCAN_RELEASE,
                      dbcan_date=DEFAULT_DBCAN_DATE, viral_loc=None, peptidase_loc=None,
                      vogdb_loc=None, vogdb_version='latest', vog_annotations=None,
                      genome_summary_form_loc=None, module_step_form_loc=None,
                      etc_module_database_loc=None, function_heatmap_form_loc=None,
                      amg_database_loc=None,
                      camper_tar_gz_loc=None,
                      skip_uniref=False, keep_database_files=False,
                      branch='master', threads=10, verbose=True):

    # TODO make this part of the plugin system
    paths_and_functions = { #just paths as of now
        "kegg_loc":{},
        "gene_ko_link_loc": {},
        "kofam_hmm_loc": {},
        "kofam_ko_list_loc": {},
        "uniref_loc": {},
        "pfam_loc": {},
        "pfam_hmm_dat": {},
        "dbcan_loc": {},
        "dbcan_fam_activities": {},
        "dbcan_subfam_ec": {},
        "vogdb_loc": {},
        "viral_loc": {},
        "peptidase_loc": {},
        "genome_summary_form_loc": {},
        "module_step_form_loc": {},
        "function_heatmap_form_loc": {},
        "amg_database_loc": {}
    }
    for i, j in locals().items():
        if i in paths_and_functions:
            paths_and_functions[i]['path'] = j
    start_time = datetime.now()
    print('%s: Database preparation started' % str(datetime.now()))

    # check inputs
    if skip_uniref and uniref_loc is not None:
        raise ValueError('Cannot skip UniRef processing and provide a location of UniRef.'
                         ' Skipping UniRef will cause provided UniRef file to not be used.')


    # check that all given files exist
    check_file_exists(*[i['path'] for i in paths_and_functions.values()])

    # setup
    if not path.isdir(output_dir):
        mkdir(output_dir)
    temporary = path.join(output_dir, 'database_files')
    mkdir(temporary)
    #setup logging
    setup_logger(os.path.join(output_dir, LOGGER))
    # output dbs dic
    output_dbs = dict()

    # Download DBs
    if vog_annotations is None:
        vog_annotations = download_vog_annotations(output_dir, vogdb_version, verbose=verbose)
    if dbcan_fam_activities is None:
        dbcan_fam_activities = download_dbcan_descriptions(
            output_dir=output_dir, dbcan_release=dbcan_version,
            upload_date=dbcan_date, verbose=verbose)
    if dbcan_subfam_ec is None:
        dbcan_subfam_ec = download_dbcan_subfam_ec(
            output_dir=output_dir, dbcan_release=dbcan_version,
            upload_date=dbcan_date, verbose=verbose)
    if dbcan_loc is None:
        dbcan_loc = download_dbcan(temporary, dbcan_release=dbcan_version, verbose=verbose)
    if pfam_hmm_dat is None:
        pfam_hmm_dat = download_pfam_descriptions(output_dir, verbose=verbose)
    if kofam_hmm_loc is None:
        kofam_hmm_loc = download_kofam_hmms(temporary, verbose=verbose)
    if genome_summary_form_loc is None:
        genome_summary_form_loc = download_genome_summary_form(temporary, branch, verbose)
    if module_step_form_loc is None:
        module_step_form_loc = download_module_step_form(temporary, branch, verbose)
    if etc_module_database_loc is None:
        etc_module_database_loc = download_etc_module_database(temporary, branch, verbose)
    if function_heatmap_form_loc is None:
        function_heatmap_form_loc = download_function_heatmap_form(temporary, branch, verbose)
    if amg_database_loc is None:
        amg_database_loc = download_amg_database(temporary, branch, verbose)
    if uniref_loc is None: 
        uniref_loc = download_uniref(temporary, uniref_version=uniref_version, verbose=verbose)
    if pfam_loc is None: 
        pfam_loc = download_pfam(temporary, verbose=verbose)
    if viral_loc is None: 
        viral_loc = download_viral_refseq(temporary, verbose=verbose)
    if peptidase_loc is None: 
        peptidase_loc = download_merops_peptidases(temporary, verbose=verbose)
    if vogdb_loc is None: 
        vogdb_loc = download_vogdb(temporary, vogdb_release=vogdb_version, verbose=verbose)
    if kofam_ko_list_loc is None: 
        kofam_ko_list_loc = download_kofam_ko_list(temporary, verbose=verbose)
    if kofam_ko_list_loc is None: 
        kofam_ko_list_loc = download_kofam_ko_list(temporary, verbose=verbose)
    if camper_tar_gz_loc is None:
        camper_tar_gz_loc = download_camper_tar_gz(temporary, verbose=verbose)

    output_dbs['genome_summary_form_loc'] = genome_summary_form_loc
    output_dbs['module_step_form_loc'] = module_step_form_loc
    output_dbs['etc_module_database_loc'] = etc_module_database_loc
    output_dbs['function_heatmap_form_loc'] = function_heatmap_form_loc
    output_dbs['amg_database_loc'] = amg_database_loc
    
    print("All raw data files where downloaded successfully")

    # Process databases

    output_dbs['dbcan_db_loc'] = process_dbcan(dbcan_loc, verbose=verbose)
    LOGGER.info('dbCAN database processed')

    if kegg_loc is not None:
        output_dbs['kegg_db_loc'] = process_kegg_db(temporary, kegg_loc, gene_ko_link_loc, kegg_download_date, threads,
                                                    verbose)
        LOGGER.info('KEGG database processed')
    if not skip_uniref:
        output_dbs['uniref_db_loc'] = process_uniref(uniref_fasta_zipped=uniref_loc, output_dir=temporary, 
                                                     uniref_version=uniref_version,
                                                     threads=threads, verbose=verbose)
        LOGGER.info('UniRef database processed')
    output_dbs['pfam_db_loc'] = download_and_process_pfam(pfam_loc, temporary,
                                                          threads=threads, verbose=verbose)
    LOGGER.info('PFAM database processed')
    output_dbs['viral_db_loc'] = download_and_process_viral_refseq(viral_loc, temporary, threads=threads,
                                                                   verbose=verbose)
    LOGGER.info('RefSeq viral database processed')
    output_dbs['peptidase_db_loc'] = download_and_process_merops_peptidases(peptidase_loc, temporary, threads=threads,
                                                                            verbose=verbose)
    LOGGER.info('MEROPS database processed')
    output_dbs['vogdb_db_loc'] = download_and_process_vogdb(vogdb_loc, temporary, vogdb_release=vogdb_version,
                                                            verbose=verbose)
    LOGGER.info('VOGdb database processed')
    output_dbs['kofam_hmm_loc'] = process_kofam_hmms(kofam_hmm_loc, temporary, verbose=verbose)
    LOGGER.info('KOfam database processed')
    output_dbs['kofam_ko_list_loc'] = download_and_process_kofam_ko_list(kofam_ko_list_loc, temporary, verbose=verbose)
    LOGGER.info('KOfam ko list processed')

    # get pfam, dbcan and vogdb descriptions
    output_dbs['pfam_hmm_dat'] = pfam_hmm_dat
    LOGGER.info('PFAM hmm dat processed')
    output_dbs['dbcan_fam_activities'] = dbcan_fam_activities
    output_dbs['dbcan_subfam_ec'] = dbcan_subfam_ec
    LOGGER.info('dbCAN fam activities processed')
    output_dbs['vog_annotations'] = vog_annotations
    LOGGER.info('VOGdb annotations processed')

    camper_locs = process_camper(camper_tar_gz_loc, 
                       temporary, threads=threads, verbose=verbose)
    output_dbs['camper_fa_db_cutoffs_loc'] = camper_locs['camper_blast_scores']
    output_dbs['camper_fa_db_loc'] = camper_locs['camper_blast']
    output_dbs['camper_hmm_loc'] = camper_locs['camper_hmm']
    output_dbs['camper_hmm_cutoffs_loc'] = camper_locs['camper_hmm_scores']
    print('%s: CAMPER annotations processed' % str(datetime.now() - start_time))

    # add genome summary form and function heatmap form
    LOGGER.info('DRAM databases and forms downloaded')

    # move all files from temporary to output that will be kept
    for db_name, output_db in output_dbs.items():
        for db_file in glob('%s*' % output_db):
            move(db_file, path.join(output_dir, path.basename(db_file)))
        output_dbs[db_name] = path.join(output_dir, path.basename(output_db))
    LOGGER.info('Files moved to final destination')

    output_dbs['description_db_loc'] = path.realpath(path.join(output_dir, 'description_db.sqlite'))

    db_handler = DatabaseHandler()
    db_handler.populate_description_db(output_dbs['description_db_loc'], update_config=False)
    db_handler.set_database_paths(**output_dbs)
    LOGGER.info('DRAM description database populated')

    if not keep_database_files:
        rmtree(temporary)
    LOGGER.info('Database preparation completed')


def update_dram_forms(output_dir, branch='master'):
    if not path.isdir(output_dir):
        mkdir(output_dir)

    form_locs = dict()
    form_locs['genome_summary_form_loc'] = download_genome_summary_form(output_dir, branch)
    form_locs['module_step_form_loc'] = download_module_step_form(output_dir, branch)
    form_locs['etc_module_database_loc'] = download_etc_module_database(output_dir, branch)
    form_locs['function_heatmap_form_loc'] = download_function_heatmap_form(output_dir, branch)
    form_locs['amg_database_loc'] = download_amg_database(output_dir, branch)
    db_handler = DatabaseHandler()
    db_handler.set_database_paths(**form_locs)



"""
import os
os.system('DRAM-setup.py -h ')
    version             print DRAM version
    prepare_databases   Download and process databases for annotation
    set_database_locations
                        Set database locations for already processed databases
    update_description_db
                        Update description database
    update_dram_forms   Update DRAM distillate and liquor forms
    print_config        Print database locations
    import_config       Import CONFIG file
    export_config       Export CONFIG file
os.system('DRAM-setup.py print_config')
os.system('DRAM-setup.py export_config')
os.system('DRAM-setup.py import_config --config_loc mag_annotator/CONFIG')
os.system('DRAM-setup.py set_database_locations')
os.sjystem('DRAM-setup.py update_description_db')
os.system('DRAM-setup.py prepare_databases --output_dir download_test --help')
os.system('DRAM-setup.py prepare_databases --output_dir download_test'
          ' --kegg_loc KEGG_LOC /home/Database/KEGG/kegg-all-orgs_20220129/kegg-all-orgs_unique_reheader_20220129.pep" ' # KEGG protein file, should be a single .pep, please merge all KEGG pep files (default: None)
          '--threads 30' Number of threads to use building mmseqs2 databases (default: 10)
          # '--gene_ko_link_loc '        # KEGG gene ko link, can be gzipped or not
          '--kofam_hmm_loc'              # hmm file for KOfam (profiles.tar.gz) (default: None)
          '--kofam_ko_list_loc'          # KOfam ko list file (ko_list.gz) (default: None)
          ' --gene_ko_link_loc'          # GENE_KO_LINK_LOC KEGG gene ko link, can be gzipped or not (default: None)
          ' --kegg_download_date'        # Date KEGG was download to include in database name (default: None)
          ' --uniref_loc'                # File path to uniref, if already downloaded (uniref90.fasta.gz) (default: None)
          ' --uniref_version'            #  UniRef version to download (default: 90)
          ' --skip_uniref'               # Do not download and process uniref90. Saves time and memory usage and does not impact DRAM distillation
          ' --pfam_loc'                  # File path to pfam-A full file, if already downloaded (Pfam-A.full.gz) (default: None)
          ' --pfam_hmm_dat'              # pfam hmm .dat file to get PF descriptions, if already downloaded (Pfam-A.hmm.dat.gz) (default: None)
          ' --dbcan_loc'                 # File path to dbCAN, if already downloaded (dbCAN-HMMdb-V9.txt) (default: None)
          ' --dbcan_fam_activities'      # CAZY family activities file, if already downloaded (CAZyDB.07302020.fam-activities.txt) (default: None)
          ' --dbcan_version'             # version of dbCAN to use (default: 10)
          ' --vogdb_loc'                 # hmm file for vogdb, if already downloaded (vog.hmm.tar.gz) (default: None)
          ' --vog_annotations'           # vogdb annotations file, if already downloaded (vog.annotations.tsv.gz) (default: None)
          ' --viral_loc'                 # File path to merged viral protein faa, if already downloaded (viral.x.protein.faa.gz) (default: None)
          ' --peptidase_loc'             # File path to MEROPS peptidase fasta, if already downloaded (pepunit.lib) (default: None)
          ' --genome_summary_form_loc'   # File path to genome summary form,if already downloaded (default: None)
          ' --module_step_form_loc'      # File path to module step form, ifalready downloaded (default: None)
          ' --etc_module_database_loc'   # File path to etc module database, if already downloaded (default: None)
          ' --function_heatmap_form_loc' # File path to function heatmap form, if already downloaded (default: None)
          ' --amg_database_loc'          # File path to amg database, if already downloaded (default: None)
          ' --branch'                    # git branch from which to download forms; THIS SHOULD NOT BE CHANGED BY REGULAR USERS (default: master)
          ' --keep_database_files'       # Keep unporcessed database files (default: False)                                                          

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
)
"""

