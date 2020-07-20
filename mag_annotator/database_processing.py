from os import path, mkdir, remove
from shutil import copy2
from datetime import datetime
from shutil import move, rmtree
from glob import glob
import json
import gzip
import tarfile
import pandas as pd
from collections import defaultdict
from skbio import read as read_sequence
from skbio import write as write_sequence

from mag_annotator.utils import run_process, make_mmseqs_db, get_database_locs, download_file, merge_files, \
    get_config_loc, remove_prefix
from mag_annotator.database_handler import DatabaseHandler
from mag_annotator.database_setup import create_description_db


# TODO: check if dbcan or pfam is down, raise appropriate error
# TODO: upgrade to pigz?


def get_iso_date():
    return datetime.today().strftime('%Y%m%d')


def check_file_exists(db_loc):
    if db_loc is None:
        return True
    elif path.isfile(db_loc):
        return True
    else:
        raise ValueError("Database location does not exist: %s" % db_loc)


def make_header_dict_from_mmseqs_db(mmseqs_db):
    mmseqs_headers_handle = open('%s_h' % mmseqs_db, 'rb')
    mmseqs_headers = mmseqs_headers_handle.read().decode(errors='ignore')
    mmseqs_headers = [i.strip() for i in mmseqs_headers.strip().split('\n\x00') if len(i) > 0]
    return [{'id': i.split(' ')[0], 'description': i} for i in mmseqs_headers]


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


def download_and_process_uniref(uniref_fasta_zipped=None, output_dir='.', uniref_version='90', threads=10,
                                verbose=True):
    """"""
    if uniref_fasta_zipped is None:  # download database if not provided
        uniref_fasta_zipped = path.join(output_dir, 'uniref%s.fasta.gz' % uniref_version)
        uniref_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz' % \
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


def process_pfam_descriptions(pfam_hmm_dat):
    if pfam_hmm_dat.endswith('.gz'):
        f = gzip.open(pfam_hmm_dat, 'r').read().decode('utf-8')
    else:
        f = open(pfam_hmm_dat).read()
    entries = f.strip().split('//')
    description_list = list()
    for i, entry in enumerate(entries):
        if len(entry) > 0:
            entry = entry.split('\n')
            ascession = None
            description = None
            for line in entry:
                line = line.strip()
                if line.startswith('#=GF AC'):
                    ascession = line.split('   ')[-1]
                if line.startswith('#=GF DE'):
                    description = line.split('   ')[-1]
            description_list.append({'id': ascession, 'description': description})
    return description_list


def download_and_process_dbcan(dbcan_hmm=None, output_dir='.', dbcan_release='8', verbose=True):
    if dbcan_hmm is None:  # download database if not provided
        dbcan_hmm = path.join(output_dir, 'dbCAN-HMMdb-V%s.txt' % dbcan_release)
        download_file('http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V%s.txt' % dbcan_release, dbcan_hmm,
                      verbose=verbose)
    run_process(['hmmpress', '-f', dbcan_hmm], verbose=verbose)
    return dbcan_hmm


def download_dbcan_descriptions(output_dir='.', upload_date='07312019', verbose=True):
    dbcan_fam_activities = path.join(output_dir, 'CAZyDB.%s.fam-activities.txt' % upload_date)
    download_file('http://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.%s.fam-activities.txt' % upload_date,
                  dbcan_fam_activities, verbose=verbose)
    return dbcan_fam_activities


def process_dbcan_descriptions(dbcan_fam_activities):
    f = open(dbcan_fam_activities)
    description_list = list()
    for line in f.readlines():
        if not line.startswith('#') and len(line.strip()) != 0:
            line = line.strip().split()
            if len(line) == 1:
                description = line[0]
            elif line[0] == line[1]:
                description = ' '.join(line[1:])
            else:
                description = ' '.join(line)
            description_list.append({'id': line[0], 'description': description.replace('\n', ' ')})
    return description_list


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


def process_vogdb_descriptions(vog_annotations):
    annotations_table = pd.read_csv(vog_annotations, sep='\t', index_col=0)
    annotations_list = [{'id': vog, 'description': '%s; %s' % (row['ConsensusFunctionalDescription'],
                                                               row['FunctionalCategory'])}
                        for vog, row in annotations_table.iterrows()]
    return annotations_list


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


def check_exists_and_add_to_location_dict(loc, name, dict_to_update):
    if loc is not None:  # if location give and exists then add to dict, else raise ValueError
        if check_file_exists(loc):
            dict_to_update[name] = path.realpath(loc)
    else:  # if location not given and is not in dict then set to none, else leave previous value
        if name not in dict_to_update:
            dict_to_update[name] = None
    return dict_to_update


def add_to_description_db(loc, name, get_description_list, db_handler):
    if loc is not None:  # if location give and exists then add to dict, else raise ValueError
        description_list = get_description_list(loc)
        db_handler.add_descriptions_to_database(description_list, name, clear_table=True)


def set_database_paths(kegg_db_loc=None, kofam_hmm_loc=None, kofam_ko_list_loc=None, uniref_db_loc=None,
                       pfam_db_loc=None, pfam_hmm_dat=None, dbcan_db_loc=None, dbcan_fam_activities=None,
                       viral_db_loc=None, peptidase_db_loc=None, vogdb_db_loc=None, vog_annotations=None,
                       description_db_loc=None, genome_summary_form_loc=None, module_step_form_loc=None,
                       etc_module_database_loc=None, function_heatmap_form_loc=None, amg_database_loc=None,
                       start_time=None, config_loc=None, use_current_locs=True, update_description_db=False):
    if start_time is None:
        start_time = datetime.now()
    print('%s: Setting database paths' % str(datetime.now() - start_time))
    if use_current_locs:
        db_dict = get_database_locs()
    else:
        db_dict = {}

    db_dict = check_exists_and_add_to_location_dict(kegg_db_loc, 'kegg', db_dict)
    db_dict = check_exists_and_add_to_location_dict(kofam_hmm_loc, 'kofam', db_dict)
    db_dict = check_exists_and_add_to_location_dict(kofam_ko_list_loc, 'kofam_ko_list', db_dict)
    db_dict = check_exists_and_add_to_location_dict(uniref_db_loc, 'uniref', db_dict)
    db_dict = check_exists_and_add_to_location_dict(pfam_db_loc, 'pfam', db_dict)
    db_dict = check_exists_and_add_to_location_dict(pfam_hmm_dat, 'pfam_hmm_dat', db_dict)
    db_dict = check_exists_and_add_to_location_dict(dbcan_db_loc, 'dbcan', db_dict)
    db_dict = check_exists_and_add_to_location_dict(dbcan_fam_activities, 'dbcan_fam_activities', db_dict)
    db_dict = check_exists_and_add_to_location_dict(viral_db_loc, 'viral', db_dict)
    db_dict = check_exists_and_add_to_location_dict(peptidase_db_loc, 'peptidase', db_dict)
    db_dict = check_exists_and_add_to_location_dict(vogdb_db_loc, 'vogdb', db_dict)
    db_dict = check_exists_and_add_to_location_dict(vog_annotations, 'vog_annotations', db_dict)

    db_dict = check_exists_and_add_to_location_dict(genome_summary_form_loc, 'genome_summary_form', db_dict)
    db_dict = check_exists_and_add_to_location_dict(module_step_form_loc, 'module_step_form', db_dict)
    db_dict = check_exists_and_add_to_location_dict(etc_module_database_loc, 'etc_module_database', db_dict)
    db_dict = check_exists_and_add_to_location_dict(function_heatmap_form_loc, 'function_heatmap_form', db_dict)
    db_dict = check_exists_and_add_to_location_dict(amg_database_loc, 'amg_database', db_dict)

    if description_db_loc is not None:
        db_dict['description_db'] = description_db_loc
    elif 'description_db' not in db_dict:
        db_dict['description_db'] = None
    print('%s: Database locations added to CONFIG' % str(datetime.now() - start_time))

    if update_description_db:
        populate_description_db(db_dict, start_time)
        check_file_exists(description_db_loc)
        print('%s: Database descriptions updated' % str(datetime.now() - start_time))

    # change data paths
    if config_loc is None:
        config_loc = get_config_loc()
    with open(config_loc, 'w') as f:
        f.write(json.dumps(db_dict))
    print('%s: Database locations set' % str(datetime.now() - start_time))


def populate_description_db(output_loc=None, db_dict=None, start_time=None):
    if start_time is None:
        start_time = datetime.now()
        print('%s: Populating description database' % str(datetime.now() - start_time))
    # setup
    if db_dict is None:
        db_dict = get_database_locs()

    if db_dict['description_db'] is None and output_loc is not None:
        db_dict['description_db'] = output_loc
    elif db_dict['description_db'] is None and output_loc is None:
        raise ValueError('Must provide output location if description db location is not set in configuration')
    elif path.exists(db_dict['description_db']):
        remove(db_dict['description_db'])

    create_description_db(db_dict['description_db'])
    db_handler = DatabaseHandler(db_dict['description_db'])
    print('%s: Database connection established' % str(datetime.now() - start_time))

    # fill database
    add_to_description_db(db_dict['kegg'], 'kegg_description', make_header_dict_from_mmseqs_db,
                          db_handler)
    print('%s: KEGG descriptions added to description database' % str(datetime.now() - start_time))
    add_to_description_db(db_dict['uniref'], 'uniref_description', make_header_dict_from_mmseqs_db,
                          db_handler)
    print('%s: UniRef descriptions added to description database' % str(datetime.now() - start_time))
    add_to_description_db(db_dict['pfam_hmm_dat'], 'pfam_description', process_pfam_descriptions,
                          db_handler)
    print('%s: PFAM descriptions added to description database' % str(datetime.now() - start_time))
    add_to_description_db(db_dict['dbcan_fam_activities'], 'dbcan_description',
                          process_dbcan_descriptions, db_handler)
    print('%s: dbCAN descriptions added to description database' % str(datetime.now() - start_time))
    add_to_description_db(db_dict['viral'], 'viral_description', make_header_dict_from_mmseqs_db,
                          db_handler)
    print('%s: RefSeq viral descriptions added to description database' % str(datetime.now() - start_time))
    add_to_description_db(db_dict['peptidase'], 'peptidase_description',
                          make_header_dict_from_mmseqs_db, db_handler)
    print('%s: MEROPS descriptions added to description database' % str(datetime.now() - start_time))
    add_to_description_db(db_dict['vog_annotations'], 'vogdb_description', process_vogdb_descriptions,
                          db_handler)
    print('%s: VOGdb descriptions added to description database' % str(datetime.now() - start_time))
    print('%s: Description database populated' % str(datetime.now() - start_time))


def prepare_databases(output_dir, kegg_loc=None, gene_ko_link_loc=None, kofam_hmm_loc=None, kofam_ko_list_loc=None,
                      kegg_download_date=None, uniref_loc=None, uniref_version='90', pfam_loc=None, pfam_hmm_dat=None,
                      dbcan_loc=None, dbcan_version='8', dbcan_fam_activities=None, dbcan_date='07312019',
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

    # get databases
    output_dbs = dict()
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
    output_dbs['dbcan_db_loc'] = download_and_process_dbcan(dbcan_loc, temporary, dbcan_release=dbcan_version,
                                                            verbose=verbose)
    print('%s: dbCAN database processed' % str(datetime.now() - start_time))
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
    if pfam_hmm_dat is None:
        pfam_hmm_dat = download_pfam_descriptions(output_dir, verbose=verbose)
    output_dbs['pfam_hmm_dat'] = pfam_hmm_dat
    print('%s: PFAM hmm dat processed' % str(datetime.now() - start_time))
    if dbcan_fam_activities is None:
        dbcan_fam_activities = download_dbcan_descriptions(output_dir, dbcan_date, verbose=verbose)
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

    set_database_paths(**output_dbs, use_current_locs=False, update_description_db=True, start_time=start_time)

    if not keep_database_files:
        rmtree(temporary)
    print('%s: Database preparation completed' % str(datetime.now() - start_time))


def print_database_locations(db_locs=None):
    if db_locs is None:
        db_locs = get_database_locs()

    print('KEGG db: %s' % db_locs.get('kegg'))
    print('KOfam db: %s' % db_locs.get('kofam'))
    print('KOfam KO list: %s' % db_locs.get('kofam_ko_list'))
    print('UniRef db: %s' % db_locs.get('uniref'))
    print('Pfam db: %s' % db_locs.get('pfam'))
    print('Pfam hmm dat: %s' % db_locs.get('pfam_hmm_dat'))
    print('dbCAN db: %s' % db_locs.get('dbcan'))
    print('dbCAN family activities: %s' % db_locs.get('dbcan_fam_activities'))
    print('RefSeq Viral db: %s' % db_locs.get('viral'))
    print('MEROPS peptidase db: %s' % db_locs.get('peptidase'))
    print('VOGDB db: %s' % db_locs.get('vogdb'))
    print('VOG annotations: %s' % db_locs.get('vog_annotations'))
    print('Description db: %s' % db_locs.get('description_db'))
    print('Genome summary form: %s' % db_locs.get('genome_summary_form'))
    print('Module step form: %s' % db_locs.get('module_step_form'))
    print('ETC module database: %s' % db_locs.get('etc_module_database'))
    print('Function heatmap form: %s' % db_locs.get('function_heatmap_form'))
    print('AMG database: %s' % db_locs.get('amg_database'))


def update_dram_forms(output_dir, branch='master'):
    if not path.isdir(output_dir):
        mkdir(output_dir)

    form_locs = dict()
    form_locs['genome_summary_form_loc'] = download_and_process_genome_summary_form(output_dir, branch)
    form_locs['module_step_form_loc'] = download_and_process_module_step_form(output_dir, branch)
    form_locs['etc_module_database_loc'] = download_and_process_etc_module_database(output_dir, branch)
    form_locs['function_heatmap_form_loc'] = download_and_process_function_heatmap_form(output_dir, branch)
    form_locs['amg_database_loc'] = download_and_process_amg_database(output_dir, branch)
    set_database_paths(**form_locs, update_description_db=False)


def export_config(output_file=None):
    if output_file is None:
        print(open(get_config_loc()).read())
    else:
        copy2(get_config_loc(), output_file)


def import_config(config_loc):
    copy2(config_loc, get_config_loc())
