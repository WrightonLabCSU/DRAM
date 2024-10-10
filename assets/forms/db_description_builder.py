"""
Contains most of the backend for the DRAM_setup.py script, used to setup databases for each user.
"""
from os import path, mkdir
from datetime import datetime
from shutil import move, rmtree, copyfile
from glob import glob
import gzip
from collections import defaultdict
import logging
import tarfile

from skbio import read as read_sequence
from skbio import write as write_sequence

from dram2.utils.database_handler import DatabaseHandler
from db_utils import (
    run_process,
    download_file,
    merge_files,
    remove_prefix,
    setup_logger,
)

NUMBER_OF_VIRAL_FILES = 2


DEFAULT_DBCAN_RELEASE = "10"
DEFAULT_DBCAN_DATE = "07292021"
DEFAULT_UNIREF_VERSION = "90"
DEFAULT_VOGDB_VERSION = "latest"
DFLT_OUTPUT_DIR = "."
LOGGER = logging.getLogger("database_processing.log")
DEFAULT_MMMSPRO_DB_NAME = "db"

from dram2.camper_kit import download as download_camper_tar_gz
from dram2.camper_kit import process as process_camper_tar_gz
from dram2.camper_kit import DOWNLOAD_OPTIONS as CAMPER_DOWNLOAD_OPTIONS
from dram2.camper_kit import PROCESS_OPTIONS as CAMPER_PROCESS_OPTIONS
from dram2.camper_kit import DRAM_SETTINGS as CAMPER_DRAM_SETTINGS
from dram2.fegenie_kit import download as download_fegenie_tar_gz
from dram2.fegenie_kit import process as process_fegenie_tar_gz
from dram2.fegenie_kit import DOWNLOAD_OPTIONS as FEGENIE_DOWNLOAD_OPTIONS
from dram2.fegenie_kit import PROCESS_OPTIONS as FEGENIE_PROCESS_OPTIONS
from dram2.fegenie_kit import DRAM_SETTINGS as FEGENIE_DRAM_SETTINGS
from dram2.sulphur_kit import download as download_sulphur_tar_gz
from dram2.sulphur_kit import process as process_sulphur_tar_gz
from dram2.sulphur_kit import DOWNLOAD_OPTIONS as SULPHUR_DOWNLOAD_OPTIONS
from dram2.sulphur_kit import PROCESS_OPTIONS as SULPHUR_PROCESS_OPTIONS
from dram2.sulphur_kit import DRAM_SETTINGS as SULPHUR_DRAM_SETTINGS

KEGG_CITATION = "Kanehisa, M., Furumichi, M., Sato, Y., Ishiguro-Watanabe, M., and Tanabe, M.; KEGG: integrating viruses and cellular organisms. Nucleic Acids Res. 49, D545-D551 (2021)."
GENE_KO_LINK_CITATION = ""
KOFAM_CITATION = ""
UNIREF_CITATION = ""
PFAM_CITATION = ""
DBCAN_CITATION = ""
VOGDB_CITATION = ""
VIRAL_REFSEQ_CITATION = ""
PEPTIDASE_CITATION = ""
DRAM_CITATION = ""
# TODO: check if dbcan or pfam is down, raise appropriate error
# TODO: upgrade to pigz?


def get_iso_date():
    return datetime.today().strftime("%Y%m%d")


def download_pfam_hmm(output_dir=".", logger=LOGGER, verbose=True):
    pfam_hmm = path.join(output_dir, "Pfam-A.hmm.dat.gz")
    link_path = (
        "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz"
    )
    logger.debug(f"Downloading Pfam from: {link_path}")
    download_file(link_path, logger, pfam_hmm, verbose=verbose)
    return pfam_hmm


def download_dbcan(
    output_dir=".",
    logger=LOGGER,
    dbcan_hmm=None,
    version=DEFAULT_DBCAN_RELEASE,
    verbose=True,
):
    dbcan_hmm = path.join(output_dir, f"dbCAN-HMMdb-V{version}.txt")
    if int(version) < int(DEFAULT_DBCAN_RELEASE):
        link_path = (
            f"http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V{version}.txt"
        )
    else:
        link_path = f"http://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V{version}.txt"

    logger.debug(f"Downloading dbCAN from: {link_path}")
    download_file(link_path, logger, dbcan_hmm, verbose=verbose)
    return dbcan_hmm


def download_dbcan_fam_activities(
    output_dir=".",
    logger=LOGGER,
    version=DEFAULT_DBCAN_RELEASE,
    upload_date=DEFAULT_DBCAN_DATE,
    verbose=True,
):
    dbcan_fam_activities = path.join(
        output_dir, f"CAZyDB.{upload_date}.fam-activities.txt"
    )
    link_path = f"https://bcb.unl.edu/dbCAN2/download/Databases/V{version}/CAZyDB.{upload_date}.fam-activities.txt"
    logger.info(f"Downloading dbCAN family activities from : {link_path}")
    download_file(link_path, logger, dbcan_fam_activities, verbose=verbose)
    return dbcan_fam_activities


def download_dbcan_subfam_ec(
    output_dir=".",
    logger=LOGGER,
    version=DEFAULT_DBCAN_RELEASE,
    upload_date=DEFAULT_DBCAN_DATE,
    verbose=True,
):
    dbcan_subfam_ec = path.join(output_dir, f"CAZyDB.{upload_date}.fam.subfam.ec.txt")
    link_path = (
        f"https://bcb.unl.edu/dbCAN2/download/Databases/"
        f"V{version}/CAZyDB.{upload_date}.fam.subfam.ec.txt"
    )
    logger.info(f"Downloading dbCAN sub-family encumber from : {link_path}")
    download_file(link_path, logger, dbcan_subfam_ec, verbose=verbose)
    return dbcan_subfam_ec


def download_kofam_hmm(output_dir=".", logger=LOGGER, verbose=False):
    kofam_profile_tar_gz = path.join(output_dir, "kofam_profiles.tar.gz")
    download_file(
        "ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz",
        logger,
        kofam_profile_tar_gz,
        verbose=verbose,
    )
    return kofam_profile_tar_gz


def generate_modified_kegg_fasta(kegg_fasta, gene_ko_link_loc=None):
    """
    Takes kegg fasta file and gene ko link file, adds kos not already in headers to headers
    Whish I knew about this, oh well I may split this out.
    """
    genes_ko_dict = defaultdict(list)
    if gene_ko_link_loc is not None:
        if gene_ko_link_loc.endswith(".gz"):
            gene_ko_link_fh = gzip.open(gene_ko_link_loc, "rt")
        else:
            gene_ko_link_fh = open(gene_ko_link_loc)
        for line in gene_ko_link_fh:
            gene, ko = line.strip().split()
            genes_ko_dict[gene].append(remove_prefix(ko, "ko:"))
    for seq in read_sequence(kegg_fasta, format="fasta"):
        new_description = seq.metadata["description"]
        for ko in genes_ko_dict[seq.metadata["id"]]:
            if ko not in new_description:
                new_description += "; %s" % ko
        seq.metadata["description"] = new_description
        yield seq


def process_kegg(
    kegg_loc,
    output_dir,
    logger,
    gene_ko_link_loc=None,
    download_date=None,
    threads=10,
    verbose=True,
):
    if download_date is None:
        download_date = get_iso_date()
    if gene_ko_link_loc is not None:
        # add KOs to end of header where KO is not already there
        kegg_mod_loc = path.join(output_dir, "kegg.mod.fa")
        write_sequence(
            generate_modified_kegg_fasta(kegg_loc, gene_ko_link_loc),
            format="fasta",
            into=kegg_mod_loc,
        )
    else:
        kegg_mod_loc = kegg_loc
    # make mmseqsdb from modified kegg fasta
    kegg_mmseqs_db = path.join(output_dir, "kegg.%s.mmsdb" % download_date)
    make_mmseqs_db(
        kegg_mod_loc,
        kegg_mmseqs_db,
        logger,
        create_index=True,
        threads=threads,
        verbose=verbose,
    )
    LOGGER.info("KEGG database processed")
    return {"kegg": kegg_mmseqs_db}


def process_kofam_hmm(
    kofam_profile_tar_gz,
    output_dir=DFLT_OUTPUT_DIR,
    logger=LOGGER,
    threads=1,
    verbose=False,
):
    kofam_profiles = path.join(output_dir, "kofam_profiles")
    mkdir(kofam_profiles)
    run_process(
        ["tar", "-xzf", kofam_profile_tar_gz, "-C", kofam_profiles],
        logger,
        verbose=verbose,
    )
    merged_kofam_profiles = path.join(output_dir, "kofam_profiles.hmm")
    merge_files(
        glob(path.join(kofam_profiles, "profiles", "*.hmm")), merged_kofam_profiles
    )
    run_process(["hmmpress", "-f", merged_kofam_profiles], logger, verbose=verbose)
    LOGGER.info("KOfam database processed")
    return {"kofam_hmm": merged_kofam_profiles}


def download_kofam_ko_list(output_dir=".", logger=LOGGER, verbose=False):
    kofam_ko_list_gz = path.join(output_dir, "kofam_ko_list.tsv.gz")
    download_file(
        "ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz",
        logger,
        kofam_ko_list_gz,
        verbose=verbose,
    )
    return kofam_ko_list_gz


def download_pfam(output_dir=".", logger=LOGGER, verbose=True):
    pfam_full_zipped = path.join(output_dir, "Pfam-A.full.gz")
    download_file(
        "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz",
        logger,
        pfam_full_zipped,
        verbose=verbose,
    )
    return pfam_full_zipped


def download_viral(
    output_dir=".", logger=LOGGER, viral_files=NUMBER_OF_VIRAL_FILES, verbose=True
):
    """Can only download newest version"""
    # download all of the viral protein files, need to know the number of files
    # TODO: Make it so that you don't need to know number of viral files in refseq viral

    faa_base_name = "viral.%s.protein.faa.gz"
    viral_faa_glob = path.join(output_dir, faa_base_name % "*")
    for number in range(viral_files):
        number += 1
        refseq_url = (
            "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.%s.protein.faa.gz"
            % number
        )
        refseq_faa = path.join(output_dir, faa_base_name % number)
        download_file(refseq_url, logger, refseq_faa, verbose=verbose)

    # then merge files from above
    merged_viral_faas = path.join(output_dir, "viral.merged.protein.faa.gz")
    run_process(
        ["cat %s > %s" % (" ".join(glob(viral_faa_glob)), merged_viral_faas)],
        logger,
        shell=True,
    )
    return merged_viral_faas


def download_uniref(
    output_dir=".",
    logger=LOGGER,
    version=DEFAULT_UNIREF_VERSION,
    threads=10,
    verbose=True,
):
    uniref_fasta_zipped = path.join(output_dir, "uniref%s.fasta.gz" % version)
    uniref_url = (
        "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz"
        % (version, version)
    )
    download_file(uniref_url, logger, uniref_fasta_zipped, verbose=verbose)
    return uniref_fasta_zipped


def download_peptidase(output_dir=".", logger=LOGGER, verbose=True):
    peptidase_faa = path.join(output_dir, "merops_peptidases_nr.faa")
    merops_url = "ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib"
    download_file(merops_url, logger, peptidase_faa, verbose=verbose)
    return peptidase_faa


def download_vogdb(
    output_dir=".", logger=LOGGER, version=DEFAULT_VOGDB_VERSION, verbose=True
):
    vog_hmm_targz = path.join(output_dir, "vog.hmm.tar.gz")
    vogdb_url = f"http://fileshare.csb.univie.ac.at/vog/{version}/vog.hmm.tar.gz"
    download_file(vogdb_url, logger, vog_hmm_targz, verbose=verbose)
    return vog_hmm_targz


def process_kofam_ko_list(
    kofam_ko_list_gz, output_dir=".", logger=LOGGER, threads=1, verbose=False
):
    # TODO: fix this so that it is gunzipped to the path
    kofam_ko_list = path.join(output_dir, "kofam_ko_list.tsv")
    run_process(
        ["gunzip", "-c", kofam_ko_list_gz],
        logger,
        save_output=kofam_ko_list,
        verbose=verbose,
    )
    LOGGER.info("KOfam ko list processed")
    return {"kofam_ko_list": kofam_ko_list}


def process_uniref(
    uniref_fasta_zipped,
    output_dir=".",
    logger=LOGGER,
    version=DEFAULT_UNIREF_VERSION,
    threads=10,
    verbose=True,
):
    uniref_mmseqs_db = path.join(
        output_dir, "uniref%s.%s.mmsdb" % (version, get_iso_date())
    )
    make_mmseqs_db(
        uniref_fasta_zipped,
        uniref_mmseqs_db,
        logger,
        create_index=True,
        threads=threads,
        verbose=verbose,
    )
    LOGGER.info("UniRef database processed")
    return {"uniref": uniref_mmseqs_db}


def process_mmspro(
    full_alignment,
    output_dir,
    logger=LOGGER,
    db_name=DEFAULT_MMMSPRO_DB_NAME,
    threads=10,
    verbose=True,
):
    mmseqs_msa = path.join(output_dir, "%s.mmsmsa" % db_name)
    run_process(
        ["mmseqs", "convertmsa", full_alignment, mmseqs_msa], logger, verbose=verbose
    )
    mmseqs_profile = path.join(output_dir, "%s.mmspro" % db_name)
    run_process(
        [
            "mmseqs",
            "msa2profile",
            mmseqs_msa,
            mmseqs_profile,
            "--match-mode",
            "1",
            "--threads",
            str(threads),
        ],
        logger,
        verbose=verbose,
    )
    tmp_dir = path.join(output_dir, "tmp")
    run_process(
        [
            "mmseqs",
            "createindex",
            mmseqs_profile,
            tmp_dir,
            "-k",
            "5",
            "-s",
            "7",
            "--threads",
            str(threads),
        ],
        logger,
        verbose=verbose,
    )
    return mmseqs_profile


def process_pfam(
    pfam_full_zipped, output_dir=".", logger=LOGGER, threads=10, verbose=True
):
    pfam_profile = process_mmspro(
        pfam_full_zipped, output_dir, logger, "pfam", threads, verbose
    )
    LOGGER.info("PFAM database processed")
    return {"pfam": pfam_profile}


def process_dbcan(input, output_dir, logger=LOGGER, verbose=True, threads=1):
    output = path.join(output_dir, path.basename(input))
    move(input, output)
    run_process(["hmmpress", "-f", output], logger, verbose=verbose)
    LOGGER.info("dbCAN database processed")
    return {"dbcan": output}


def process_viral(
    merged_viral_faas,
    output_dir=".",
    logger=LOGGER,
    viral_files=NUMBER_OF_VIRAL_FILES,
    threads=10,
    verbose=True,
):
    refseq_viral_mmseqs_db = path.join(
        output_dir, "refseq_viral.%s.mmsdb" % get_iso_date()
    )
    make_mmseqs_db(
        merged_viral_faas,
        refseq_viral_mmseqs_db,
        logger,
        create_index=True,
        threads=threads,
        verbose=verbose,
    )
    LOGGER.info("RefSeq viral database processed")
    return {"viral": refseq_viral_mmseqs_db}


def process_peptidase(
    peptidase_faa, output_dir=".", logger=LOGGER, threads=10, verbose=True
):
    peptidase_mmseqs_db = path.join(output_dir, "peptidases.%s.mmsdb" % get_iso_date())
    make_mmseqs_db(
        peptidase_faa,
        peptidase_mmseqs_db,
        logger,
        create_index=True,
        threads=threads,
        verbose=verbose,
    )
    LOGGER.info("MEROPS database processed")
    return {"peptidase": peptidase_mmseqs_db}


def process_vogdb(
    vog_hmm_targz,
    output_dir=".",
    logger=LOGGER,
    version=DEFAULT_VOGDB_VERSION,
    threads=1,
    verbose=True,
):
    hmm_dir = path.join(output_dir, "vogdb_hmms")
    mkdir(hmm_dir)
    vogdb_targz = tarfile.open(vog_hmm_targz)
    vogdb_targz.extractall(hmm_dir)
    vog_hmms = path.join(output_dir, f"vog_{version}_hmms.txt")
    merge_files(glob(path.join(hmm_dir, "VOG*.hmm")), vog_hmms)
    run_process(["hmmpress", "-f", vog_hmms], logger, verbose=verbose)
    LOGGER.info("VOGdb database processed")
    return {"vogdb": vog_hmms}


def download_vog_annotations(
    output_dir, logger=LOGGER, version=DEFAULT_VOGDB_VERSION, verbose=True
):
    vog_annotations = path.join(output_dir, "vog_annotations_%s.tsv.gz" % version)
    download_file(
        "http://fileshare.csb.univie.ac.at/vog/%s/vog.annotations.tsv.gz" % version,
        logger,
        vog_annotations,
        verbose=verbose,
    )
    return vog_annotations


def download_genome_summary_form(
    output_dir, logger=LOGGER, branch="master", verbose=True
):
    genome_summary_form = path.join(
        output_dir, "genome_summary_form.%s.tsv" % get_iso_date()
    )
    download_file(
        "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/genome_summary_form.tsv"
        % branch,
        logger,
        genome_summary_form,
        verbose=verbose,
    )
    return genome_summary_form


def download_module_step_form(output_dir, logger=LOGGER, branch="master", verbose=True):
    function_heatmap_form = path.join(
        output_dir, "module_step_form.%s.tsv" % get_iso_date()
    )
    download_file(
        "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/module_step_form.tsv"
        % branch,
        logger,
        function_heatmap_form,
        verbose=verbose,
    )
    return function_heatmap_form


def download_etc_module_database(
    output_dir, logger=LOGGER, branch="master", verbose=True
):
    etc_module_database = path.join(
        output_dir, "etc_mdoule_database.%s.tsv" % get_iso_date()
    )
    download_file(
        "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/etc_module_database.tsv"
        % branch,
        logger,
        etc_module_database,
        verbose=verbose,
    )
    return etc_module_database


def download_function_heatmap_form(
    output_dir, logger=LOGGER, branch="master", verbose=True
):
    function_heatmap_form = path.join(
        output_dir, "function_heatmap_form.%s.tsv" % get_iso_date()
    )
    download_file(
        "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/function_heatmap_form.tsv"
        % branch,
        logger,
        function_heatmap_form,
        verbose=verbose,
    )
    return function_heatmap_form


def download_amg_database(output_dir, logger=LOGGER, branch="master", verbose=True):
    amg_database = path.join(output_dir, "amg_database.%s.tsv" % get_iso_date())
    download_file(
        "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/amg_database.tsv"
        % branch,
        logger,
        amg_database,
        verbose=verbose,
    )
    return amg_database


def check_file_exists(*paths):
    for i in paths:
        if i is None:
            continue
        elif path.isfile(i):
            continue
        else:
            raise ValueError(f"Database location does not exist: {i}")


def prepare_databases(
    output_dir,
    loggpath=None,
    kegg_loc=None,
    gene_ko_link_loc=None,
    kofam_hmm_loc=None,
    kofam_ko_list_loc=None,
    kegg_download_date=None,
    uniref_loc=None,
    uniref_version=DEFAULT_UNIREF_VERSION,
    pfam_loc=None,
    pfam_hmm_loc=None,
    dbcan_loc=None,
    dbcan_fam_activities: str = None,
    dbcan_subfam_ec: str = None,
    dbcan_version=DEFAULT_DBCAN_RELEASE,
    dbcan_date=DEFAULT_DBCAN_DATE,
    viral_loc=None,
    peptidase_loc=None,
    vogdb_loc=None,
    vogdb_version=DEFAULT_VOGDB_VERSION,
    vog_annotations=None,
    genome_summary_form_loc=None,
    module_step_form_loc=None,
    etc_module_database_loc=None,
    function_heatmap_form_loc=None,
    amg_database_loc=None,
    camper_tar_gz_loc=None,
    number_of_viral_files=NUMBER_OF_VIRAL_FILES,
    fegenie_tar_gz_loc=None,
    sulphur_tar_gz_loc=None,
    skip_uniref=False,
    keep_database_files=False,
    branch="master",
    threads=10,
    verbose=True,
    select_db=None,
    clear_config=False,
):

    dram_settings = {
        "kegg": {
            "name": "KEGG db",
            "description_db_updated": "Unknown, or Never",
            "citation": KEGG_CITATION,
        },
        "gene_ko_link": {"name": "KEGG Gene KO link", "citation": KEGG_CITATION},
        "kofam_hmm": {"name": "KOfam db", "citation": KOFAM_CITATION},
        "kofam_ko_list": {"name": "KOfam KO list", "citation": KOFAM_CITATION},
        "uniref": {
            "name": "UniRef db",
            "description_db_updated": "Unknown, or Never",
            "citation": UNIREF_CITATION,
        },
        "pfam": {"name": "Pfam db", "citation": PFAM_CITATION},
        "pfam_hmm": {
            "name": "Pfam hmm dat",
            "description_db_updated": "Unknown, or Never",
            "citation": PFAM_CITATION,
        },
        "dbcan": {"name": "dbCAN db", "citation": DBCAN_CITATION},
        "dbcan_fam_activities": {
            "name": "dbCAN family activities",
            "citation": DBCAN_CITATION,
        },
        "dbcan_subfam_ec": {
            "name": "dbCAN subfamily EC numbers",
            "citation": DBCAN_CITATION,
        },
        "vogdb": {"name": "VOGDB db", "citation": VOGDB_CITATION},
        "vog_annotations": {
            "name": "VOG annotations",
            "description_db_updated": "Unknown, or Never",
            "citation": VOGDB_CITATION,
        },
        "viral": {
            "name": "RefSeq Viral db",
            "description_db_updated": "Unknown, or Never",
            "citation": VIRAL_REFSEQ_CITATION,
        },
        "peptidase": {
            "name": "MEROPS peptidase db",
            "description_db_updated": "Unknown, or Never",
            "citation": PEPTIDASE_CITATION,
        },
        "genome_summary_form": {"name": "Genome summary form"},
        "module_step_form": {"name": "Module step form"},
        "function_heatmap_form": {"name": "Function heatmap form"},
        "amg_database": {"name": "AMG database"},
        "etc_module_database": {"name": "ETC module database"},
    }
    dram_settings.update(CAMPER_DRAM_SETTINGS)
    dram_settings.update(FEGENIE_DRAM_SETTINGS)
    dram_settings.update(SULPHUR_DRAM_SETTINGS)
    database_settings = {
        "kegg": {},
        "gene_ko_link": {},
        "kofam_hmm": {},
        "kofam_ko_list": {},
        "uniref": {"version": uniref_version},
        "pfam": {},
        "pfam_hmm": {},
        "dbcan": {"version": dbcan_version},
        "dbcan_fam_activities": {"version": dbcan_version, "upload_date": dbcan_date},
        "dbcan_subfam_ec": {"version": dbcan_version, "upload_date": dbcan_date},
        "vogdb": {"version": vogdb_version},
        "vog_annotations": {"version": vogdb_version},
        "viral": {},
        "peptidase": {},
        "genome_summary_form": {"branch": branch},
        "module_step_form": {"branch": branch},
        "function_heatmap_form": {"branch": branch},
        "amg_database": {"branch": branch},
        "etc_module_database": {"branch": branch},
    }
    database_settings.update(CAMPER_DOWNLOAD_OPTIONS)
    database_settings.update(FEGENIE_DOWNLOAD_OPTIONS)
    database_settings.update(SULPHUR_DOWNLOAD_OPTIONS)
    process_settings = {
        "kegg": {},
        "gene_ko_link": {},
        "kofam_hmm": {},
        "kofam_ko_list": {},
        "uniref": {"version": uniref_version},
        "pfam": {},
        "pfam_hmm": {},
        "dbcan": {},
        "dbcan_fam_activities": {},
        "dbcan_subfam_ec": {},
        "vogdb": {},
        "vog_annotations": {},
        "viral": {"viral_files": number_of_viral_files},
        "peptidase": {},
        "genome_summary_form": {},
        "module_step_form": {},
        "function_heatmap_form": {},
        "amg_database": {},
        "etc_module_database": {},
    }
    process_settings.update(CAMPER_PROCESS_OPTIONS)
    process_settings.update(FEGENIE_PROCESS_OPTIONS)
    process_settings.update(SULPHUR_PROCESS_OPTIONS)

    # setup temp, logging, and db_handler
    if not path.isdir(output_dir):
        mkdir(output_dir)
    temporary = path.join(output_dir, "database_files")
    mkdir(temporary)
    main_log = path.join(output_dir, "database_processing.log")
    setup_logger(LOGGER, *[(main_log, loggpath) if loggpath is not None else main_log])
    db_handler = DatabaseHandler(logger=LOGGER)
    if clear_config or select_db is None:
        db_handler.clear_config()

    db_handler.config["log_path"] = main_log
    db_handler.write_config()
    LOGGER.info("Starting the process of downloading data")

    if skip_uniref:
        LOGGER.info("Skipping UniRef")
        del database_settings["uniref"]

    locs = {
        i.removesuffix("_loc"): j
        for i, j in locals().items()
        if i.endswith("_loc") and j is not None
    }
    download_functions = {
        i.removeprefix("download_"): j
        for i, j in globals().items()
        if callable(j) and i.startswith("download_")
    }
    process_functions = {
        i.removeprefix("process_"): j
        for i, j in globals().items()
        if callable(j) and i.startswith("process_")
    }
    functions = {
        i: j for i, j in globals().items() if callable(j) and i.startswith("download_")
    }

    # Check any specified paths exist
    missing_user_inputs = [i for i in locs if not path.exists(i)]
    if len(missing_user_inputs) > 1:
        raise ValueError(
            f"The fallowing user provided paths don't seem to exist: {missing_user_inputs}"
        )

    un_obtainable = [
        i for i in database_settings if i not in locs and i not in download_functions
    ]

    for i in un_obtainable:
        LOGGER.info(
            f"The {i}_loc argument was not used to specify a downloaded {i} file, and dram can not"
            " download it its self. So it is assumed that the user wants to set up DRAM without it"
        )
        del database_settings[i]

    # check inputs
    if skip_uniref and uniref_loc is not None:
        raise ValueError(
            "Cannot skip UniRef processing and provide a location of UniRef."
            " Skipping UniRef will cause provided UniRef file to not be used."
        )

    if select_db is not None:
        miss_name = [i for i in select_db if i not in database_settings]
        user_inputs = [i for i in locs if i not in select_db]
        if len(miss_name) > 0:
            LOGGER.error(
                "Only the databases in the db list can be pased to select_db, "
                f"you passed {miss_name} which is/are not in the list."
            )
            raise ValueError("Bad user input, see log")

        if len(user_inputs) > 0:
            LOGGER.error(
                f"The user provided location for {user_inputs}, but required it not be used by proving"
                f" the select_db argument for other databases. This would suggest that the"
                " user may have made a mistake and so this error is rased."
            )
            raise ValueError("Bad user input, see log")

        database_settings = {i: database_settings[i] for i in select_db}

    LOGGER.info("Database preparation started")

    # Download DBs
    for i, j in database_settings.items():
        if locs.get(i) is None:
            LOGGER.info(f"Downloading {i}")
            if i in process_functions:
                locs[i] = download_functions[i](temporary, LOGGER, **j, verbose=verbose)
            else:
                locs[i] = download_functions[i](
                    output_dir, LOGGER, **j, verbose=verbose
                )
            j["Download time"] = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
            j["Origin"] = "Downloaded by DRAM"
        else:
            j = {k: "Unknown" for k in j}
            j["Download time"] = "Unknown"
            j["Origin"] = "Provided by user"
            j["Original path"] = locs[i]
            if i not in process_functions:
                LOGGER.info(f"Copying {locs[i]} to output_dir")
                locs[i] = copyfile(
                    locs[i], path.join(output_dir, path.basename(locs[i]))
                )

    LOGGER.info("All raw data files were downloaded successfully")

    # Process databases
    for i in locs:
        processed_locs = {}
        if i in process_functions:
            LOGGER.info(f"Processing {i}")
            processed_locs = process_functions[i](
                locs[i],
                output_dir,
                LOGGER,
                threads=threads,
                verbose=verbose,
                **process_settings[i],
            )
        else:
            processed_locs = {i: locs[i]}
        for k, v in processed_locs.items():
            final_dest = path.join(output_dir, path.basename(v))
            if v != final_dest:
                for db_file in glob("%s*" % v):
                    move(db_file, path.join(output_dir, path.basename(db_file)))
                v = path.join(output_dir, path.basename(v))
            # update_dram_forms the settings per OUTPUT fill, including the process_settings
            #  and database_settings, which are per input file.
            if db_handler.config.get("setup_info") is None:
                db_handler.config["setup_info"] = {}
            db_handler.config["setup_info"][k] = {
                **dram_settings[k],
                **process_settings[i],
                **database_settings[i],
            }
            db_handler.set_database_paths(**{f"{k}_loc": v})
            db_handler.write_config()
            LOGGER.info(f"Moved {k} to final destination, configuration updated")

    LOGGER.info("Populating the description db, this may take some time")
    db_handler.config["description_db"] = path.realpath(
        path.join(output_dir, "description_db.sqlite")
    )
    db_handler.write_config()
    db_handler.populate_description_db(
        db_handler.config["description_db"], select_db, update_config=False
    )
    # todo make db handler such that the destruction on success write_config
    db_handler.write_config()
    LOGGER.info("DRAM description database populated")

    if not keep_database_files:
        rmtree(temporary)
    LOGGER.info("Database preparation completed")


def update_dram_forms(output_dir, branch="master"):
    if not path.isdir(output_dir):
        mkdir(output_dir)
    form_locs = dict()
    form_locs["genome_summary_form_loc"] = download_genome_summary_form(
        output_dir, branch
    )
    form_locs["module_step_form_loc"] = download_module_step_form(output_dir, branch)
    form_locs["etc_module_database_loc"] = download_etc_module_database(
        output_dir, branch
    )
    form_locs["function_heatmap_form_loc"] = download_function_heatmap_form(
        output_dir, branch
    )
    form_locs["amg_database_loc"] = download_amg_database(output_dir, branch)
    db_handler = DatabaseHandler()
    db_handler.set_database_paths(**form_locs)


"""

os.system("DRAM.py annotate -i /home/projects-wrighton-2/DRAM/development_flynn/release_validation/data_sets/15_soil_genomes/all_data/Cytophaga_hutchinsonii_ATCC_33406.fasta  -o test_15soil --use_camper --use_fegenie")
import os
os.system('DRAM-setup.py prepare_databases --output_dir /home/projects-wrighton-2/DRAM/dram_data/dram1.4_final_06_21_22/ --kegg_loc /home/Database/KEGG/kegg-all-orgs_20220603/kegg-all-orgs_unique_reheader_20220603.pep --threads 40')
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
os.system('DRAM-setup.py import_config --config_loc dram2/CONFIG')
os.system('DRAM-setup.py set_database_locations')
os.system('DRAM-setup.py update_description_db')

os.system('DRAM-setup.py prepare_databases --output_dir download_test')
os.system('DRAM-setup.py prepare_databases --output_dir download_test'
os.system('DRAM-setup.py prepare_databases --output_dir download_test --select_db  vogdb')
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
          ' --pfam_hmm'              # pfam hmm .dat file to get PF descriptions, if already downloaded (Pfam-A.hmm.dat.gz) (default: None)
          ' --dbcan_loc'                 # File path to dbCAN, if already downloaded (dbCAN-HMMdb-V9.txt) (default: None)
          ' --dbcan_fam_activities'      # CAZY family activities file, if already downloaded (CAZyDB.07302020.fam-activities.txt) (default: None)
          ' --dbcan_version'             # version of dbCAN to use (default: 10)
          ' --vogdb_loc'                 # hmm file for vogdb, if already downloaded (vog.hmm.tar.gz) (default: None)
          ' --vog_annotations_loc'       # vogdb annotations file, if already downloaded (vog.annotations.tsv.gz) (default: None)
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