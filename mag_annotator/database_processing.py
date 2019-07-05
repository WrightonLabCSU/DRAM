from os import path, mkdir
from datetime import datetime
from shutil import move, rmtree
from glob import glob
from pkg_resources import resource_filename
import json
import gzip

from mag_annotator.utils import run_process, make_mmseqs_db, get_database_locs, download_file


def get_iso_date():
    return datetime.today().strftime('%Y%m%d')


def check_file_exists(db_loc):
    if path.isfile(db_loc):
        return True
    else:
        raise ValueError("Database location does not exist: %s" % db_loc)


def process_kegg_db(output_dir, kegg_loc, download_date=None, threads=10, verbose=True):
    check_file_exists(kegg_loc)
    if download_date is None:
        download_date = get_iso_date()
    kegg_mmseqs_db = path.join(output_dir, 'kegg.%s.mmsdb' % download_date)
    make_mmseqs_db(kegg_loc, kegg_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
    return kegg_mmseqs_db


def download_and_process_unifref(uniref_fasta_zipped=None, output_dir='.', uniref_version='90', threads=10,
                                 verbose=True):
    """"""
    if uniref_fasta_zipped is None:  # download database if not provided
        uniref_fasta_zipped = path.join(output_dir, 'uniref%s.fasta.gz' % uniref_version)
        uniref_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz' %\
                     (uniref_version, uniref_version)
        download_file(uniref_url, uniref_fasta_zipped, verbose=verbose)
    else:
        check_file_exists(uniref_fasta_zipped)
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
    else:
        check_file_exists(pfam_full_zipped)
    pfam_profile = process_mmspro(pfam_full_zipped, output_dir, 'pfam', threads, verbose)
    return pfam_profile


def download_pfam_descriptions(output_dir='.', pfam_release='32.0', verbose=True):
    pfam_hmm_dat = path.join(output_dir, 'Pfam-A.hmm.dat.gz')
    download_file('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam%s/Pfam-A.hmm.dat.gz' % pfam_release,
                  pfam_hmm_dat, verbose=verbose)
    return pfam_hmm_dat


def process_pfam_descriptions(pfam_hmm_dat=None):
    check_file_exists(pfam_hmm_dat)
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
    else:
        check_file_exists(dbcan_hmm)
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
    else:
        check_file_exists(merged_viral_faas)

    # make mmseqs database
    refseq_viral_mmseqs_db = path.join(output_dir, 'refseq_viral.%s.mmsdb' % get_iso_date())
    make_mmseqs_db(merged_viral_faas, refseq_viral_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
    return refseq_viral_mmseqs_db


def download_and_process_merops_peptidases(peptidase_faa=None, output_dir='.', threads=10, verbose=True):
    if peptidase_faa is None:  # download database if not provided
        peptidase_faa = path.join(output_dir, 'merops_peptidases_nr.faa')
        merops_url = 'ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib'
        download_file(merops_url, peptidase_faa, verbose=verbose)
    else:
        check_file_exists(peptidase_faa)
    peptidase_mmseqs_db = path.join(output_dir, 'peptidases.%s.mmsdb' % get_iso_date())
    make_mmseqs_db(peptidase_faa, peptidase_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
    return peptidase_mmseqs_db


def download_and_process_kegg_modules(output_dir):
    module_summary_form_path = path.join(output_dir, 'module_summary_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/shafferm/checkMetab/master/data/module_summary_form.tsv',
                  module_summary_form_path, verbose=True)
    return module_summary_form_path


def download_and_process_genome_summary_form(output_dir):
    genome_summary_form = path.join(output_dir, 'genome_summary_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/shafferm/checkMetab/master/data/genome_summary_form.tsv',
                  genome_summary_form, verbose=True)
    return genome_summary_form


def check_exists_and_add_to_location_dict(loc, name, dict_to_update):
    if loc is not None:  # if location give and exists then add to dict, else raise ValueError
        if check_file_exists(loc):
            dict_to_update[name] = path.abspath(loc)
    else:  # if location not given and is not in dict then set to none, else leave previous value
        if name not in dict_to_update:
            dict_to_update[name] = None
    return dict_to_update


def set_database_paths(kegg_db_loc=None, uniref_db_loc=None, pfam_db_loc=None, pfam_hmm_dat=None, dbcan_db_loc=None,
                       viral_db_loc=None, peptidase_db_loc=None,
                       module_step_form_loc=None, genome_summary_form_loc=None):
    """Processes pfam_hmm_dat"""
    db_dict = get_database_locs()
    db_dict = check_exists_and_add_to_location_dict(kegg_db_loc, 'kegg', db_dict)
    db_dict = check_exists_and_add_to_location_dict(uniref_db_loc, 'uniref', db_dict)
    db_dict = check_exists_and_add_to_location_dict(pfam_db_loc, 'pfam', db_dict)
    if pfam_hmm_dat is None and 'pfam_description' not in db_dict:
        db_dict['pfam_description'] = None
    elif check_file_exists(pfam_hmm_dat):
        db_dict['pfam_description'] = process_pfam_descriptions(pfam_hmm_dat)
    else:
        pass
    db_dict = check_exists_and_add_to_location_dict(dbcan_db_loc, 'dbcan', db_dict)
    db_dict = check_exists_and_add_to_location_dict(viral_db_loc, 'viral', db_dict)
    db_dict = check_exists_and_add_to_location_dict(peptidase_db_loc, 'peptidase', db_dict)
    db_dict = check_exists_and_add_to_location_dict(module_step_form_loc, 'module_step_form', db_dict)
    db_dict = check_exists_and_add_to_location_dict(genome_summary_form_loc, 'genome_summary_form', db_dict)

    # change data paths
    with open(path.abspath(resource_filename('mag_annotator', 'CONFIG')), 'w') as f:
        f.write(json.dumps(db_dict))


def prepare_databases(output_dir, kegg_loc=None, kegg_download_date=None, uniref_loc=None, uniref_version='90',
                      pfam_loc=None, pfam_release='32.0', pfam_hmm_dat=None, dbcan_loc=None, dbcan_version='7',
                      viral_loc=None, peptidase_loc=None, keep_database_files=False, threads=10, verbose=True):
    mkdir(output_dir)
    temporary = path.join(output_dir, 'database_files')
    mkdir(temporary)
    output_dbs = dict()
    if kegg_loc is not None:
        output_dbs['kegg_db_loc'] = process_kegg_db(temporary, kegg_loc, kegg_download_date, threads, verbose)
    output_dbs['uniref_db_loc'] = download_and_process_unifref(uniref_loc, temporary, uniref_version=uniref_version,
                                                               threads=threads, verbose=verbose)
    output_dbs['pfam_db_loc'] = download_and_process_pfam(pfam_loc, temporary, pfam_release=pfam_release,
                                                          threads=threads, verbose=verbose)
    output_dbs['dbcan_db_loc'] = download_and_process_dbcan(dbcan_loc, temporary, dbcan_release=dbcan_version,
                                                            verbose=verbose)
    output_dbs['viral_db_loc'] = download_and_process_viral_refseq(viral_loc, temporary, threads=threads,
                                                                   verbose=verbose)
    output_dbs['peptidase_db_loc'] = download_and_process_merops_peptidases(peptidase_loc, temporary, threads=threads,
                                                                            verbose=verbose)
    # get module step form
    output_dbs['module_step_form_loc'] = download_and_process_kegg_modules(temporary)

    # add genome summary form
    output_dbs['genome_summary_form_loc'] = download_and_process_genome_summary_form(temporary)

    for db_name, output_db in output_dbs.items():
        for db_file in glob('%s*' % output_db):
            move(db_file, path.join(output_dir, path.basename(db_file)))
        output_dbs[db_name] = path.join(output_dir, path.basename(output_db))

    # get pfam descriptions
    if pfam_hmm_dat is None:
        pfam_hmm_dat = download_pfam_descriptions(output_dir, pfam_release=pfam_release, verbose=verbose)
    output_dbs['pfam_hmm_dat'] = pfam_hmm_dat

    set_database_paths(**output_dbs)

    if not keep_database_files:
        rmtree(temporary)


def is_db_in_dict(key, dict_):
    if key in dict_:
        return dict_[key]
    else:
        return str(None)


def print_database_locations():
    db_locs = get_database_locs()

    print('KEGG db loc: %s' % is_db_in_dict('kegg', db_locs))
    print('UniRef db loc: %s' % is_db_in_dict('uniref', db_locs))
    print('Pfam db loc: %s' % is_db_in_dict('pfam', db_locs))
    has_pfam_desc = 'pfam_description' in db_locs
    print('Has Pfam descriptions: %s' % has_pfam_desc)
    print('dbCAN db loc: %s' % is_db_in_dict('dbcan', db_locs))
    print('RefSeq Viral db loc: %s' % is_db_in_dict('viral', db_locs))
    print('MEROPS peptidase db loc: %s' % is_db_in_dict('peptidase', db_locs))
    has_module_step_form = 'module_step_form' in db_locs
    print('Has module steps form: %s' % has_module_step_form)
    has_genome_summary_form = 'genome_summary_form' in db_locs
    print('Has genome summary form: %s' % has_genome_summary_form)
