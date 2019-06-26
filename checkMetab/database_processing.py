from os import path, mkdir
from checkMetab.utils import run_process, merge_files, make_mmseqs_db
from datetime import datetime
from shutil import move, rmtree
from glob import glob
from pkg_resources import resource_filename
import json
import gzip


def get_iso_date():
    return datetime.today().strftime('%Y%m%d')


def download_file(url, output_file, verbose=True):
    if verbose:
        print('downloading %s' % url)
    run_process(['wget', '-O', output_file, url], verbose=verbose)


def download_and_process_unifref(uniref_fasta_zipped=None, output_dir='.', uniref_version='90', threads=10,
                                 verbose=True):
    """"""
    if uniref_fasta_zipped is None:  # download database if not provided
        uniref_fasta_zipped = path.join(output_dir, 'uniref%s.fasta.gz' % uniref_version)
        uniref_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz' %\
                     (uniref_version, uniref_version)
        download_file(uniref_url, uniref_fasta_zipped, verbose=verbose)
    uniref_mmseqs_db = path.join(output_dir, 'uniref%s.%s.mmsdb' % (uniref_version, get_iso_date()))
    make_mmseqs_db(uniref_fasta_zipped, uniref_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
    return uniref_mmseqs_db


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


def download_and_process_pfam(pfam_full_zipped=None, output_dir='.', pfam_release='32.0', threads=10, verbose=True):
    if pfam_full_zipped is None:  # download database if not provided
        pfam_full_zipped = path.join(output_dir, 'Pfam-A.full.gz')
        download_file('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam%s/Pfam-A.full.gz' % pfam_release,
                      pfam_full_zipped)
    pfam_profile = process_mmspro(pfam_full_zipped, output_dir, 'pfam', threads, verbose)
    return pfam_profile


def download_and_process_pfam_descriptions(pfam_hmm_dat=None, output_dir='.', pfam_release='32.0', verbose=True):
    if pfam_hmm_dat is None:
        pfam_hmm_dat = path.join(output_dir, 'Pfam-A.hmm.dat.gz')
        download_file('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam%s/Pfam-A.hmm.dat.gz' % pfam_release,
                      pfam_hmm_dat, verbose=verbose)
    if pfam_hmm_dat.endswith('.gz'):
        f = gzip.open(pfam_hmm_dat, 'r').read().decode('utf-8')
    else:
        f = open(pfam_hmm_dat).read()
    entries = f.split('//')
    description_dict = dict()
    for entry in entries:
        entry = entry.split('\n')
        ascession = None
        description = None
        for line in entry:
            line = line.strip()
            if line.startswith('#=GF AC'):
                ascession = line.split('   ')[-1]
            if line.startswith('#=GF DE'):
                description = line.split('   ')[-1]
        description_dict[ascession] = description
    return description_dict


def download_and_process_dbcan(dbcan_hmm=None, output_dir='.', dbcan_release='7', verbose=True):
    if dbcan_hmm is None:  # download database if not provided
        dbcan_hmm = path.join(output_dir, 'dbCAN-HMMdb-V%s.txt' % dbcan_release)
        download_file('http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V%s.txt' % dbcan_release, dbcan_hmm,
                      verbose=verbose)
    run_process(['hmmpress', '-f', dbcan_hmm], verbose=verbose)
    return dbcan_hmm


def download_and_process_viral_refseq(merged_viral_faas=None, output_dir='.', viral_files=3, threads=10, verbose=True):
    """Can only download newest version"""
    # download all of the viral protein files, need to know the number of files
    # TODO: Make it so that you don't need to know number of viral files in refseq viral

    if merged_viral_faas is None:  # download database if not provided
        faa_base_name = 'viral.%s.protein.faa.gz'
        viral_faa_glob = path.join(output_dir, faa_base_name % '*')
        for number in range(viral_files):
            number += 1
            refseq_url = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.%s.protein.faa.gz' % number
            refseq_faa = path.join(output_dir, faa_base_name % number)
            download_file(refseq_url, refseq_faa)

        # then merge files from above
        merged_viral_faas = path.join(output_dir, 'viral.merged.protein.faa.gz')
        run_process(['cat %s > %s' % (' '.join(glob(viral_faa_glob)), merged_viral_faas)], shell=True)

    # make mmseqs database
    refseq_viral_mmseqs_db = path.join(output_dir, 'refseq_viral.%s.mmsdb' % get_iso_date())
    make_mmseqs_db(merged_viral_faas, refseq_viral_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
    return refseq_viral_mmseqs_db


def process_kegg_db(output_dir, kegg_loc, download_date=None, threads=10, verbose=True):
    if download_date is None:
        download_date = get_iso_date()
    kegg_mmseqs_db = path.join(output_dir, 'kegg.%s.mmsdb' % download_date)
    make_mmseqs_db(kegg_loc, kegg_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
    return kegg_mmseqs_db


def update_config(output_dbs):
    # change data paths
    with open(path.abspath(resource_filename('checkMetab', 'DATABASE_LOCATIONS')), 'w') as f:
        f.write(json.dumps(output_dbs))


def prepare_databases(output_dir, kegg_loc=None, kegg_download_date=None, uniref_loc=None, uniref_version=90,
                      pfam_loc=None, pfam_version=32.0, pfam_hmm_dat=None, dbcan_loc=None, dbcan_version=7,
                      viral_loc=None, keep_database_files=False, threads=10, verbose=True):
    mkdir(output_dir)
    temporary = path.join(output_dir, 'database_files')
    mkdir(temporary)
    output_dbs = dict()
    if kegg_loc is not None:
        output_dbs['kegg'] = process_kegg_db(temporary, kegg_loc, kegg_download_date, threads, verbose)
    output_dbs['uniref'] = download_and_process_unifref(uniref_loc, temporary, threads=threads, verbose=verbose)
    output_dbs['pfam'] = download_and_process_pfam(pfam_loc, temporary, threads=threads, verbose=verbose)
    output_dbs['pfam_description'] = download_and_process_pfam_descriptions(pfam_hmm_dat)
    output_dbs['dbcan'] = download_and_process_dbcan(dbcan_loc, temporary, verbose=verbose)
    output_dbs['viral'] = download_and_process_viral_refseq(viral_loc, temporary, threads=threads, verbose=verbose)

    for output_db in output_dbs:
        for db_file in glob('%s*' % output_db):
            move(db_file, path.join(output_dir, path.basename(db_file)))

    if not keep_database_files:
        rmtree(temporary)


def check_file_exists(db_loc):
    if path.isfile(db_loc):
        return True
    else:
        raise ValueError("Database location does not exist: %s" % db_loc)


def set_database_paths(kegg_db_loc=None, uniref_db_loc=None, pfam_db_loc=None, pfam_hmm_dat=None, dbcan_db_loc=None,
                       viral_db_loc=None):
    db_dict = dict()
    if kegg_db_loc is not None:
        if check_file_exists(kegg_db_loc):
            db_dict['kegg'] = kegg_db_loc
    if uniref_db_loc is not None:
        if check_file_exists(uniref_db_loc):
            db_dict['uniref'] = uniref_db_loc
    if pfam_db_loc is not None:
        if check_file_exists(pfam_db_loc):
            db_dict['pfam'] = pfam_db_loc
    if pfam_hmm_dat is not None:
        if check_file_exists(pfam_hmm_dat):
            db_dict['pfam_description'] = download_and_process_pfam_descriptions(pfam_hmm_dat)
    if dbcan_db_loc is not None:
        if check_file_exists(dbcan_db_loc):
            db_dict['dbcan'] = dbcan_db_loc
    if viral_db_loc is not None:
        if check_file_exists(viral_db_loc):
            db_dict['viral'] = viral_db_loc

    update_config(db_dict)
