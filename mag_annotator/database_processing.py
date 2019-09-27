from os import path, mkdir, remove
from datetime import datetime
from shutil import move, rmtree
from glob import glob
from pkg_resources import resource_filename
import json
import gzip
import tarfile
import pandas as pd
from collections import defaultdict
from skbio import read as read_sequence
from skbio import write as write_sequence

from mag_annotator.utils import run_process, make_mmseqs_db, get_database_locs, download_file, merge_files
from mag_annotator.database_handler import DatabaseHandler
from mag_annotator.database_setup import create_description_db

# TODO: check if dbcan or pfam is down, raise appropriate error


def get_iso_date():
    return datetime.today().strftime('%Y%m%d')


def check_file_exists(db_loc):
    if path.isfile(db_loc):
        return True
    else:
        raise ValueError("Database location does not exist: %s" % db_loc)


def make_header_dict_from_mmseqs_db(mmseqs_db):
    mmseqs_headers_handle = open('%s_h' % mmseqs_db, 'rb')
    mmseqs_headers = mmseqs_headers_handle.read().decode(errors='ignore')
    mmseqs_headers = [i.strip() for i in mmseqs_headers.strip().split(' \n\x00') if len(i) > 0]
    return [{'id': i.split(' ')[0], 'description': i} for i in mmseqs_headers]


def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text  # or whatever


def generate_modified_kegg_fasta(kegg_fasta, gene_ko_link_loc):
    """Takes kegg fasta file and gene ko link file, adds kos not already in headers to headers"""
    if gene_ko_link_loc.endswith('.gz'):
        gene_ko_link_fh = gzip.open(gene_ko_link_loc, 'rt')
    else:
        gene_ko_link_fh = open(gene_ko_link_loc)
    genes_ko_dict = defaultdict(list)
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
    check_file_exists(kegg_loc)
    check_file_exists(gene_ko_link_loc)
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


def download_and_process_uniref(uniref_fasta_zipped=None, output_dir='.', uniref_version='90', threads=10,
                                verbose=True):
    """"""
    if uniref_fasta_zipped is None:  # download database if not provided
        uniref_fasta_zipped = path.join(output_dir, 'uniref%s.fasta.gz' % uniref_version)
        uniref_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz' % \
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


def process_pfam_descriptions(pfam_hmm_dat):
    check_file_exists(pfam_hmm_dat)
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


def download_and_process_dbcan(dbcan_hmm=None, output_dir='.', dbcan_release='7', verbose=True):
    if dbcan_hmm is None:  # download database if not provided
        dbcan_hmm = path.join(output_dir, 'dbCAN-HMMdb-V%s.txt' % dbcan_release)
        download_file('http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V%s.txt' % dbcan_release, dbcan_hmm,
                      verbose=verbose)
    else:
        check_file_exists(dbcan_hmm)
    run_process(['hmmpress', '-f', dbcan_hmm], verbose=verbose)
    return dbcan_hmm


def download_dbcan_descriptions(output_dir='.', upload_date='07312018', verbose=True):
    dbcan_fam_activities = path.join(output_dir, 'CAZyDB.%s.fam-activities.txt' % upload_date)
    download_file('http://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.%s.fam-activities.txt' % upload_date,
                  dbcan_fam_activities, verbose=verbose)
    return dbcan_fam_activities


def process_dbcan_descriptions(dbcan_fam_activities):
    check_file_exists(dbcan_fam_activities)
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


def download_and_process_vogdb(vog_hmm_targz=None, output_dir='.', vogdb_release='latest', verbose=True):
    if vog_hmm_targz is None:
        vog_hmm_targz = path.join(output_dir, 'vog.hmm.tar.gz')
        vogdb_url = 'http://fileshare.csb.univie.ac.at/vog/%s/vog.hmm.tar.gz' % vogdb_release
        download_file(vogdb_url, vog_hmm_targz, verbose=verbose)
    else:
        check_file_exists(vog_hmm_targz)
    hmm_dir = path.join(output_dir, 'vogdb_hmms')
    mkdir(hmm_dir)
    vogdb_targz = tarfile.open(vog_hmm_targz)
    vogdb_targz.extractall(hmm_dir)
    vog_hmms = path.join(output_dir, 'vog_%s_hmms.txt' % vogdb_release)
    merge_files(path.join(hmm_dir, 'VOG*.hmm'), vog_hmms)
    run_process(['hmmpress', '-f', vog_hmms], verbose=verbose)
    return vog_hmms


def download_vog_annotations(output_dir, vogdb_version='latest', verbose=True):
    vog_annotations = path.join(output_dir, 'vog_annotations_%s.tsv.gz' % vogdb_version)
    download_file('http://fileshare.csb.univie.ac.at/vog/%s/vog.annotations.tsv.gz' % vogdb_version,
                  vog_annotations, verbose=verbose)
    return vog_annotations


def process_vogdb_descriptions(vog_annotations):
    check_file_exists(vog_annotations)
    annotations_table = pd.read_csv(vog_annotations, sep='\t', index_col=0)
    annotations_list = [{'id': vog, 'description': '%s; %s' % (row['ConsensusFunctionalDescription'],
                                                               row['FunctionalCategory'])}
                        for vog, row in annotations_table.iterrows()]
    return annotations_list


def download_and_process_genome_summary_form(output_dir, branch='master'):
    genome_summary_form = path.join(output_dir, 'genome_summary_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/shafferm/DRAM/%s/data/genome_summary_form.tsv' % branch,
                  genome_summary_form, verbose=True)
    return genome_summary_form


def download_and_process_module_step_form(output_dir, branch='master'):
    function_heatmap_form = path.join(output_dir, 'module_step_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/shafferm/DRAM/%s/data/module_step_form.tsv' % branch,
                  function_heatmap_form, verbose=True)
    return function_heatmap_form


def download_and_process_function_heatmap_form(output_dir, branch='master'):
    function_heatmap_form = path.join(output_dir, 'function_heatmap_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/shafferm/DRAM/%s/data/function_heatmap_form.tsv' % branch,
                  function_heatmap_form, verbose=True)
    return function_heatmap_form


def download_and_process_amg_database(output_dir, branch='master'):
    function_heatmap_form = path.join(output_dir, 'function_heatmap_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/shafferm/DRAM/%s/data/function_heatmap_form.tsv' % branch,
                  function_heatmap_form, verbose=True)
    return function_heatmap_form


def check_exists_and_add_to_location_dict(loc, name, dict_to_update):
    if loc is not None:  # if location give and exists then add to dict, else raise ValueError
        if check_file_exists(loc):
            dict_to_update[name] = path.abspath(loc)
    else:  # if location not given and is not in dict then set to none, else leave previous value
        if name not in dict_to_update:
            dict_to_update[name] = None
    return dict_to_update


def check_exists_and_add_to_description_db(loc, name, get_description_list, db_handler):
    if loc is not None:  # if location give and exists then add to dict, else raise ValueError
        if check_file_exists(loc):
            description_list = get_description_list(loc)
            db_handler.add_descriptions_to_database(description_list, name, clear_table=True)


def set_database_paths(kegg_db_loc=None, uniref_db_loc=None, pfam_db_loc=None, pfam_hmm_dat=None, dbcan_db_loc=None,
                       dbcan_fam_activities=None, viral_db_loc=None, peptidase_db_loc=None, vogdb_db_loc=None,
                       vog_annotations=None, description_db_loc=None, genome_summary_form_loc=None,
                       module_step_form_loc=None, function_heatmap_form_loc=None, amg_database_loc=None,
                       update_description_db=False):
    """Processes pfam_hmm_dat"""
    db_dict = get_database_locs()

    db_dict = check_exists_and_add_to_location_dict(kegg_db_loc, 'kegg', db_dict)
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
    db_dict = check_exists_and_add_to_location_dict(function_heatmap_form_loc, 'function_heatmap_form', db_dict)
    db_dict = check_exists_and_add_to_location_dict(amg_database_loc, 'amg_database', db_dict)

    if description_db_loc is not None:
        db_dict['description_db'] = description_db_loc

    if update_description_db:
        populate_description_db(db_dict)

    # change data paths
    with open(path.abspath(resource_filename('mag_annotator', 'CONFIG')), 'w') as f:
        f.write(json.dumps(db_dict))


def populate_description_db(db_dict=None):
    # setup
    if db_dict is None:
        db_dict = get_database_locs()
    if path.exists(db_dict['description_db']):
        remove(db_dict['description_db'])
    create_description_db(db_dict['description_db'])
    db_handler = DatabaseHandler(db_dict['description_db'])

    # fill database
    check_exists_and_add_to_description_db(db_dict['kegg'], 'kegg_description', make_header_dict_from_mmseqs_db,
                                           db_handler)
    check_exists_and_add_to_description_db(db_dict['uniref'], 'uniref_description', make_header_dict_from_mmseqs_db,
                                           db_handler)
    check_exists_and_add_to_description_db(db_dict['pfam_hmm_dat'], 'pfam_description', process_pfam_descriptions,
                                           db_handler)
    check_exists_and_add_to_description_db(db_dict['dbcan_fam_activities'], 'dbcan_description',
                                           process_dbcan_descriptions, db_handler)
    check_exists_and_add_to_description_db(db_dict['viral'], 'viral_description', make_header_dict_from_mmseqs_db,
                                           db_handler)
    check_exists_and_add_to_description_db(db_dict['peptidase'], 'peptidase_description',
                                           make_header_dict_from_mmseqs_db, db_handler)
    check_exists_and_add_to_description_db(db_dict['vog_annotations'], 'vogdb_description', process_vogdb_descriptions,
                                           db_handler)


def prepare_databases(output_dir, kegg_loc=None, gene_ko_link_loc=None, kegg_download_date=None, uniref_loc=None,
                      uniref_version='90', pfam_loc=None, pfam_release='32.0', pfam_hmm_dat=None, dbcan_loc=None,
                      dbcan_version='7', dbcan_fam_activities=None, dbcan_date='07312018', viral_loc=None,
                      peptidase_loc=None, vogdb_loc=None, vogdb_version='latest', vog_annotations=None,
                      genome_summary_form_loc=None, module_step_form_loc=None, function_heatmap_form_loc=None,
                      amg_database_loc=None, keep_database_files=False, branch='master', threads=10, verbose=True):
    # check that all given files exist
    if kegg_loc is not None:
        check_file_exists(kegg_loc)
    if gene_ko_link_loc is not None:
        check_file_exists(gene_ko_link_loc)
    if uniref_loc is not None:
        check_file_exists(uniref_loc)
    if pfam_loc is not None:
        check_file_exists(pfam_loc)
    if pfam_hmm_dat is not None:
        check_file_exists(pfam_hmm_dat)
    if dbcan_loc is not None:
        check_file_exists(dbcan_loc)
    if dbcan_fam_activities is not None:
        check_file_exists(dbcan_fam_activities)
    if vogdb_loc is not None:
        check_file_exists(vogdb_loc)
    if viral_loc is not None:
        check_file_exists(viral_loc)
    if peptidase_loc is not None:
        check_file_exists(peptidase_loc)
    if genome_summary_form_loc is not None:
        check_file_exists(genome_summary_form_loc)
    if module_step_form_loc is not None:
        check_file_exists(module_step_form_loc)
    if function_heatmap_form_loc is not None:
        check_file_exists(function_heatmap_form_loc)
    if amg_database_loc is not None:
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
    output_dbs['uniref_db_loc'] = download_and_process_uniref(uniref_loc, temporary, uniref_version=uniref_version,
                                                              threads=threads, verbose=verbose)
    output_dbs['pfam_db_loc'] = download_and_process_pfam(pfam_loc, temporary, pfam_release=pfam_release,
                                                          threads=threads, verbose=verbose)
    output_dbs['dbcan_db_loc'] = download_and_process_dbcan(dbcan_loc, temporary, dbcan_release=dbcan_version,
                                                            verbose=verbose)
    output_dbs['viral_db_loc'] = download_and_process_viral_refseq(viral_loc, temporary, threads=threads,
                                                                   verbose=verbose)
    output_dbs['peptidase_db_loc'] = download_and_process_merops_peptidases(peptidase_loc, temporary, threads=threads,
                                                                            verbose=verbose)
    output_dbs['vogdb_db_loc'] = download_and_process_vogdb(vogdb_loc, temporary, vogdb_release=vogdb_version,
                                                            verbose=verbose)

    # add genome summary form and function heatmap form
    if genome_summary_form_loc is None:
        output_dbs['genome_summary_form_loc'] = download_and_process_genome_summary_form(temporary, branch)
    else:
        output_dbs['genome_summary_form_loc'] = genome_summary_form_loc
    if module_step_form_loc is None:
        output_dbs['module_step_form_loc'] = download_and_process_module_step_form(temporary, branch)
    else:
        output_dbs['module_step_form_loc'] = module_step_form_loc
    if function_heatmap_form_loc is None:
        output_dbs['function_heatmap_form_loc'] = download_and_process_function_heatmap_form(temporary, branch)
    else:
        output_dbs['function_heatmap_form_loc'] = function_heatmap_form_loc
    if function_heatmap_form_loc is None:
        output_dbs['amg_database_loc'] = download_and_process_amg_database(temporary, branch)
    else:
        output_dbs['amg_database_loc'] = amg_database_loc

    for db_name, output_db in output_dbs.items():
        for db_file in glob('%s*' % output_db):
            move(db_file, path.join(output_dir, path.basename(db_file)))
        output_dbs[db_name] = path.join(output_dir, path.basename(output_db))

    # get pfam and dbcan descriptions
    if pfam_hmm_dat is None:
        pfam_hmm_dat = download_pfam_descriptions(output_dir, pfam_release=pfam_release, verbose=verbose)
    output_dbs['pfam_hmm_dat'] = pfam_hmm_dat
    if dbcan_fam_activities is None:
        dbcan_fam_activities = download_dbcan_descriptions(output_dir, dbcan_date, verbose=verbose)
    output_dbs['dbcan_fam_activities'] = dbcan_fam_activities
    if vog_annotations is None:
        vog_annotations = download_vog_annotations(output_dir, vogdb_version, verbose=verbose)
    output_dbs['vog_annotations'] = vog_annotations

    output_dbs['description_db_loc'] = path.abspath(path.join(output_dir, 'description_db.sqlite'))

    set_database_paths(**output_dbs, update_description_db=True)

    if not keep_database_files:
        rmtree(temporary)


def is_db_in_dict(key, dict_):
    if key in dict_:
        return dict_[key]
    else:
        return str(None)


def print_database_locations():
    db_locs = get_database_locs()

    print('KEGG db location: %s' % is_db_in_dict('kegg', db_locs))
    print('UniRef db location: %s' % is_db_in_dict('uniref', db_locs))
    print('Pfam db location: %s' % is_db_in_dict('pfam', db_locs))
    print('dbCAN db location: %s' % is_db_in_dict('dbcan', db_locs))
    print('RefSeq Viral db location: %s' % is_db_in_dict('viral', db_locs))
    print('MEROPS peptidase db location: %s' % is_db_in_dict('peptidase', db_locs))
    print('VOGDB db location: %s' % is_db_in_dict('vogdb', db_locs))
    print('Description db location: %s' % is_db_in_dict('description_db', db_locs))
    print('Genome summary form location: %s' % is_db_in_dict('genome_summary_form', db_locs))
    print('Function heatmap form location: %s' % is_db_in_dict('function_heatmap_form', db_locs))
    print('AMG database location: %s' % is_db_in_dict('amg_database', db_locs))
