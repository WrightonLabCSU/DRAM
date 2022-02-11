from os import path, mkdir
from datetime import datetime
from shutil import move, rmtree
from glob import glob
import gzip
import tarfile
from collections import defaultdict
from skbio import read as read_sequence
from skbio import write as write_sequence

from mag_annotator.database_handler import DatabaseHandler
from mag_annotator.utils import run_process, make_mmseqs_db, download_file, merge_files, remove_prefix

DEFAULT_DBCAN_RELEASE = '10'
DEFAULT_DBCAN_DATE = '07292021'
DEFAULT_UNIREF_VERSION = '90'

# TODO: check if dbcan or pfam is down, raise appropriate error
# TODO: upgrade to pigz?


def get_iso_date():
    return datetime.today().strftime('%Y%m%d')


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


def download_and_process_kofam_hmms(kofam_profile_tar_gz=None, output_dir='.', verbose=False):
    if kofam_profile_tar_gz is None:
        kofam_profile_tar_gz = path.join(output_dir, 'kofam_profiles.tar.gz')
        download_file('ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz', kofam_profile_tar_gz, verbose=verbose)
    kofam_profiles = path.join(output_dir, 'kofam_profiles')
    mkdir(kofam_profiles)
    run_process(['tar', '-xzf', kofam_profile_tar_gz, '-C', kofam_profiles], verbose=verbose)
    merged_kofam_profiles = path.join(output_dir, 'kofam_profiles.hmm')
    merge_files(glob(path.join(kofam_profiles, 'profiles', '*.hmm')), merged_kofam_profiles)
    run_process(['hmmpress', '-f', merged_kofam_profiles], verbose=verbose)
    return merged_kofam_profiles


def download_and_process_kofam_ko_list(kofam_ko_list_gz=None, output_dir='.', verbose=False):
    if kofam_ko_list_gz is None:
        kofam_ko_list_gz = path.join(output_dir, 'kofam_ko_list.tsv.gz')
        download_file('ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz', kofam_ko_list_gz, verbose=verbose)
    # TODO: fix this so that it is gunzipped to the path
    kofam_ko_list = path.join(output_dir, 'kofam_ko_list.tsv')
    run_process(['gunzip', kofam_ko_list_gz], verbose=verbose)
    return kofam_ko_list


def download_and_process_uniref(uniref_fasta_zipped=None, output_dir='.', uniref_version=DEFAULT_UNIREF_VERSION, threads=10,
                                verbose=True):
    """"""
    if uniref_fasta_zipped is None:  # download database if not provided
        uniref_fasta_zipped = path.join(output_dir, 'uniref%s.fasta.gz' % uniref_version)
        uniref_url = 'https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz' % \
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


def download_and_process_pfam(pfam_full_zipped=None, output_dir='.', threads=10, verbose=True):
    if pfam_full_zipped is None:  # download database if not provided
        pfam_full_zipped = path.join(output_dir, 'Pfam-A.full.gz')
        download_file('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz', pfam_full_zipped,
                      verbose=verbose)
    pfam_profile = process_mmspro(pfam_full_zipped, output_dir, 'pfam', threads, verbose)
    return pfam_profile


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


def process_dbcan(dbcan_hmm, verbose=True):
    run_process(['hmmpress', '-f', dbcan_hmm], verbose=verbose)
    return dbcan_hmm


def download_dbcan_descriptions(output_dir='.', dbcan_release=DEFAULT_DBCAN_RELEASE, upload_date=DEFAULT_DBCAN_DATE, verbose=True):
    dbcan_fam_activities = path.join(output_dir, 'CAZyDB.%s.fam-activities.txt' % upload_date)
    link_path = f"https://bcb.unl.edu/dbCAN2/download/Databases/V{dbcan_release}/CAZyDB.{upload_date}.fam-activities.txt"
    print(f"Downloading dbCAN family activities from : {link_path}")
    download_file(link_path, dbcan_fam_activities, verbose=verbose)
    return dbcan_fam_activities


def download_and_process_viral_refseq(merged_viral_faas=None, output_dir='.', viral_files=2, threads=10, verbose=True):
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
            download_file(refseq_url, refseq_faa, verbose=verbose)

        # then merge files from above
        merged_viral_faas = path.join(output_dir, 'viral.merged.protein.faa.gz')
        run_process(['cat %s > %s' % (' '.join(glob(viral_faa_glob)), merged_viral_faas)], shell=True)

    # make mmseqs database
    refseq_viral_mmseqs_db = path.join(output_dir, 'refseq_viral.%s.mmsdb' % get_iso_date())
    make_mmseqs_db(merged_viral_faas, refseq_viral_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
    return refseq_viral_mmseqs_db


def download_and_process_merops_peptidases(peptidase_faa=None, output_dir='.', threads=10, verbose=True):
    if peptidase_faa is None:  # download database if not provided
        peptidase_faa = path.join(output_dir, 'merops_peptidases_nr.faa')
        merops_url = 'ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib'
        download_file(merops_url, peptidase_faa, verbose=verbose)
    peptidase_mmseqs_db = path.join(output_dir, 'peptidases.%s.mmsdb' % get_iso_date())
    make_mmseqs_db(peptidase_faa, peptidase_mmseqs_db, create_index=True, threads=threads, verbose=verbose)
    return peptidase_mmseqs_db


def download_and_process_vogdb(vog_hmm_targz=None, output_dir='.', vogdb_release='latest', verbose=True):
    if vog_hmm_targz is None:
        vog_hmm_targz = path.join(output_dir, 'vog.hmm.tar.gz')
        vogdb_url = 'http://fileshare.csb.univie.ac.at/vog/%s/vog.hmm.tar.gz' % vogdb_release
        download_file(vogdb_url, vog_hmm_targz, verbose=verbose)
    hmm_dir = path.join(output_dir, 'vogdb_hmms')
    mkdir(hmm_dir)
    vogdb_targz = tarfile.open(vog_hmm_targz)
    vogdb_targz.extractall(hmm_dir)
    vog_hmms = path.join(output_dir, 'vog_%s_hmms.txt' % vogdb_release)
    merge_files(glob(path.join(hmm_dir, 'VOG*.hmm')), vog_hmms)
    run_process(['hmmpress', '-f', vog_hmms], verbose=verbose)
    return vog_hmms


def download_vog_annotations(output_dir, vogdb_version='latest', verbose=True):
    vog_annotations = path.join(output_dir, 'vog_annotations_%s.tsv.gz' % vogdb_version)
    download_file('http://fileshare.csb.univie.ac.at/vog/%s/vog.annotations.tsv.gz' % vogdb_version,
                  vog_annotations, verbose=verbose)
    return vog_annotations


def download_and_process_genome_summary_form(output_dir, branch='master', verbose=True):
    genome_summary_form = path.join(output_dir, 'genome_summary_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/shafferm/DRAM/%s/data/genome_summary_form.tsv' % branch,
                  genome_summary_form, verbose=verbose)
    return genome_summary_form


def download_and_process_module_step_form(output_dir, branch='master', verbose=True):
    function_heatmap_form = path.join(output_dir, 'module_step_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/shafferm/DRAM/%s/data/module_step_form.tsv' % branch,
                  function_heatmap_form, verbose=verbose)
    return function_heatmap_form


def download_and_process_etc_module_database(output_dir, branch='master', verbose=True):
    etc_module_database = path.join(output_dir, 'etc_mdoule_database.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/shafferm/DRAM/%s/data/etc_module_database.tsv' % branch,
                  etc_module_database, verbose=verbose)
    return etc_module_database


def download_and_process_function_heatmap_form(output_dir, branch='master', verbose=True):
    function_heatmap_form = path.join(output_dir, 'function_heatmap_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/shafferm/DRAM/%s/data/function_heatmap_form.tsv' % branch,
                  function_heatmap_form, verbose=verbose)
    return function_heatmap_form


def download_and_process_amg_database(output_dir, branch='master', verbose=True):
    amg_database = path.join(output_dir, 'amg_database.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/shafferm/DRAM/%s/data/amg_database.tsv' % branch,
                  amg_database, verbose=verbose)
    return amg_database


def check_file_exists(db_loc):
    if db_loc is None:
        return True
    elif path.isfile(db_loc):
        return True
    else:
        raise ValueError("Database location does not exist: %s" % db_loc)


def prepare_databases(output_dir, kegg_loc=None, gene_ko_link_loc=None, kofam_hmm_loc=None, kofam_ko_list_loc=None,
                      kegg_download_date=None, uniref_loc=None, uniref_version=DEFAULT_UNIREF_VERSION, pfam_loc=None, pfam_hmm_dat=None,
                      dbcan_loc=None, dbcan_version=DEFAULT_DBCAN_RELEASE, dbcan_fam_activities=None, dbcan_date=DEFAULT_DBCAN_DATE,
                      viral_loc=None, peptidase_loc=None, vogdb_loc=None, vogdb_version='latest', vog_annotations=None,
                      genome_summary_form_loc=None, module_step_form_loc=None, etc_module_database_loc=None,
                      function_heatmap_form_loc=None, amg_database_loc=None, skip_uniref=False,
                      keep_database_files=False, branch='master', threads=10, verbose=True):
    start_time = datetime.now()
    print('%s: Database preparation started' % str(datetime.now()))

    # check inputs
    if skip_uniref and uniref_loc is not None:
        raise ValueError('Cannot skip UniRef processing and provide a location of UniRef. Skipping UniRef will cause '
                         'provided UniRef file to not be used.')

    # check that all given files exist
    check_file_exists(kegg_loc)
    check_file_exists(gene_ko_link_loc)
    check_file_exists(kofam_hmm_loc)
    check_file_exists(kofam_ko_list_loc)
    check_file_exists(uniref_loc)
    check_file_exists(pfam_loc)
    check_file_exists(pfam_hmm_dat)
    check_file_exists(dbcan_loc)
    check_file_exists(dbcan_fam_activities)
    check_file_exists(vogdb_loc)
    check_file_exists(viral_loc)
    check_file_exists(peptidase_loc)
    check_file_exists(genome_summary_form_loc)
    check_file_exists(module_step_form_loc)
    check_file_exists(function_heatmap_form_loc)
    check_file_exists(amg_database_loc)

    # setup
    if not path.isdir(output_dir):
        mkdir(output_dir)
    temporary = path.join(output_dir, 'database_files')
    mkdir(temporary)

    # Download DBs
    if dbcan_fam_activities is None:
        dbcan_fam_activities = download_dbcan_descriptions(
            output_dir=output_dir, dbcan_release=dbcan_version,
            upload_date=dbcan_date, verbose=verbose)
    if pfam_hmm_dat is None:
        pfam_hmm_dat = download_pfam_descriptions(output_dir, verbose=verbose)
    if dbcan_loc is None:
        dbcan_loc = download_dbcan(temporary, dbcan_release=dbcan_version, verbose=verbose)

    # Process databases
    output_dbs = dict()

    output_dbs['dbcan_db_loc'] = process_dbcan(dbcan_loc, verbose=verbose)
    print('%s: dbCAN database processed' % str(datetime.now() - start_time))

    if kegg_loc is not None:
        output_dbs['kegg_db_loc'] = process_kegg_db(temporary, kegg_loc, gene_ko_link_loc, kegg_download_date, threads,
                                                    verbose)
        print('%s: KEGG database processed' % str(datetime.now() - start_time))
    if not skip_uniref:
        output_dbs['uniref_db_loc'] = download_and_process_uniref(uniref_loc, temporary, uniref_version=uniref_version,
                                                                  threads=threads, verbose=verbose)
        print('%s: UniRef database processed' % str(datetime.now() - start_time))
    output_dbs['pfam_db_loc'] = download_and_process_pfam(pfam_loc, temporary,
                                                          threads=threads, verbose=verbose)
    print('%s: PFAM database processed' % str(datetime.now() - start_time))
    output_dbs['viral_db_loc'] = download_and_process_viral_refseq(viral_loc, temporary, threads=threads,
                                                                   verbose=verbose)
    print('%s: RefSeq viral database processed' % str(datetime.now() - start_time))
    output_dbs['peptidase_db_loc'] = download_and_process_merops_peptidases(peptidase_loc, temporary, threads=threads,
                                                                            verbose=verbose)
    print('%s: MEROPS database processed' % str(datetime.now() - start_time))
    output_dbs['vogdb_db_loc'] = download_and_process_vogdb(vogdb_loc, temporary, vogdb_release=vogdb_version,
                                                            verbose=verbose)
    print('%s: VOGdb database processed' % str(datetime.now() - start_time))
    output_dbs['kofam_hmm_loc'] = download_and_process_kofam_hmms(kofam_hmm_loc, temporary, verbose=verbose)
    print('%s: KOfam database processed' % str(datetime.now() - start_time))
    output_dbs['kofam_ko_list_loc'] = download_and_process_kofam_ko_list(kofam_ko_list_loc, temporary, verbose=verbose)
    print('%s: KOfam ko list processed' % str(datetime.now() - start_time))

    # get pfam, dbcan and vogdb descriptions
    output_dbs['pfam_hmm_dat'] = pfam_hmm_dat
    print('%s: PFAM hmm dat processed' % str(datetime.now() - start_time))
    output_dbs['dbcan_fam_activities'] = dbcan_fam_activities
    print('%s: dbCAN fam activities processed' % str(datetime.now() - start_time))
    if vog_annotations is None:
        vog_annotations = download_vog_annotations(output_dir, vogdb_version, verbose=verbose)
    output_dbs['vog_annotations'] = vog_annotations
    print('%s: VOGdb annotations processed' % str(datetime.now() - start_time))

    # add genome summary form and function heatmap form
    if genome_summary_form_loc is None:
        output_dbs['genome_summary_form_loc'] = download_and_process_genome_summary_form(temporary, branch, verbose)
    else:
        output_dbs['genome_summary_form_loc'] = genome_summary_form_loc
    if module_step_form_loc is None:
        output_dbs['module_step_form_loc'] = download_and_process_module_step_form(temporary, branch, verbose)
    else:
        output_dbs['module_step_form_loc'] = module_step_form_loc
    if etc_module_database_loc is None:
        output_dbs['etc_module_database_loc'] = download_and_process_etc_module_database(temporary, branch, verbose)
    else:
        output_dbs['etc_module_database_loc'] = etc_module_database_loc
    if function_heatmap_form_loc is None:
        output_dbs['function_heatmap_form_loc'] = download_and_process_function_heatmap_form(temporary, branch, verbose)
    else:
        output_dbs['function_heatmap_form_loc'] = function_heatmap_form_loc
    if amg_database_loc is None:
        output_dbs['amg_database_loc'] = download_and_process_amg_database(temporary, branch, verbose)
    else:
        output_dbs['amg_database_loc'] = amg_database_loc
    print('%s: DRAM databases and forms downloaded' % str(datetime.now() - start_time))

    # move all files from temporary to output that will be kept
    for db_name, output_db in output_dbs.items():
        for db_file in glob('%s*' % output_db):
            move(db_file, path.join(output_dir, path.basename(db_file)))
        output_dbs[db_name] = path.join(output_dir, path.basename(output_db))
    print('%s: Files moved to final destination' % str(datetime.now() - start_time))

    output_dbs['description_db_loc'] = path.realpath(path.join(output_dir, 'description_db.sqlite'))

    db_handler = DatabaseHandler()
    db_handler.populate_description_db(output_dbs['description_db_loc'], update_config=False)
    db_handler.set_database_paths(**output_dbs)
    print('%s: DRAM description database populated' % str(datetime.now() - start_time))

    if not keep_database_files:
        rmtree(temporary)
    print('%s: Database preparation completed' % str(datetime.now() - start_time))


def update_dram_forms(output_dir, branch='master'):
    if not path.isdir(output_dir):
        mkdir(output_dir)

    form_locs = dict()
    form_locs['genome_summary_form_loc'] = download_and_process_genome_summary_form(output_dir, branch)
    form_locs['module_step_form_loc'] = download_and_process_module_step_form(output_dir, branch)
    form_locs['etc_module_database_loc'] = download_and_process_etc_module_database(output_dir, branch)
    form_locs['function_heatmap_form_loc'] = download_and_process_function_heatmap_form(output_dir, branch)
    form_locs['amg_database_loc'] = download_and_process_amg_database(output_dir, branch)
    db_handler = DatabaseHandler()
    db_handler.set_database_paths(**form_locs)
