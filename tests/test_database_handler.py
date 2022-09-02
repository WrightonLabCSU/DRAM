import pytest

from os import path, environ
from sqlalchemy.orm.exc import NoResultFound
import json

from mag_annotator.database_setup import KeggDescription, create_description_db
from mag_annotator.database_handler import DatabaseHandler, print_database_locations, import_config, export_config, get_config_loc
from mag_annotator.utils import setup_logger
import logging
import pandas as pd
import numpy as np

TEST_CONF_PATH = path.join('tests', 'data', 'test_CONFIG')
TEST_CONF_PATH_1_3 = path.join('tests', 'data', 'test_CONFIG_1.3-')


@pytest.fixture()
def logger(tmpdir):
    logger = logging.getLogger('test_log')
    setup_logger(logger)
    return logger


@pytest.fixture()
def db_handler(logger):
    return DatabaseHandler(logger, TEST_CONF_PATH)


def test_db_hander_init():  # TODO: this
    pass

def test_add_descriptions_to_database(tmpdir, logger):
    db_handler = DatabaseHandler(logger, path.join('tests', 'data', 'test_CONFIG'))
    test_db = path.join(tmpdir.mkdir('test_db'), 'test_db.sqlite')
    create_description_db(test_db)
    db_handler.config['description_db'] = test_db
    db_handler.start_db_session()
    kegg_entries = [{'id': 'K00001', 'description': 'The first KO'},
                    {'id': 'K00002', 'description': 'The second KO'}]
    db_handler.add_descriptions_to_database(kegg_entries, 'kegg_description')
    assert len(db_handler.session.query(KeggDescription).all()) == 2
    # now test that clear table works
    kegg_entry = [{'id': 'K00003', 'description': 'The third KO'}]
    db_handler.add_descriptions_to_database(kegg_entry, 'kegg_description', clear_table=True)
    assert len(db_handler.session.query(KeggDescription).all()) == 1


@pytest.fixture()
def db_w_entries(tmpdir, logger):
    db_handler = DatabaseHandler(logger, path.join('tests', 'data', 'test_CONFIG'))
    test_db = path.join(tmpdir.mkdir('test_db_w_entries'), 'test_db.sqlite')
    create_description_db(test_db)
    db_handler.config['description_db'] = test_db
    db_handler.start_db_session()
    db_handler.add_descriptions_to_database([
        {'id': 'K00001', 'description': 'The first KO'},
        {'id': 'K00002', 'description': 'The second KO'}], 
                                            'kegg_description',
                                            clear_table=True)
    return db_handler


def test_get_description(db_w_entries):
    description = db_w_entries.get_description('K00001', 'kegg_description')
    assert description == 'The first KO'
    with pytest.raises(NoResultFound):
        db_w_entries.get_description('K00003', 'kegg_description')
    description_dict = db_w_entries.get_descriptions(['K00001', 'K00002'], 'kegg_description')
    assert isinstance(description_dict, dict)
    assert len(description_dict) == 2
    description_dict = db_w_entries.get_descriptions(['K00001', 'K00002', 'K00003'], 'kegg_description')
    assert isinstance(description_dict, dict)
    assert len(description_dict) == 2
    description_dict = db_w_entries.get_descriptions(['K00003'], 'kegg_description')
    assert isinstance(description_dict, dict)
    assert len(description_dict) == 0


def test_get_database_names(db_w_entries):
    names = db_w_entries.get_database_names()
    assert len(names) == 7


def test_filter_db_locs(logger):
    db_handler = DatabaseHandler(logger, path.join('tests', 'data', 'test_CONFIG'))
    db_handler.config['search_databases'] = {
        'kegg': 'a/fake/location/kegg.fasta', 
        'kofam_hmm': 'a/fake/location/kofam.hmm',
        'kofam_ko_list': 'a/fake/location/kofam_ko_list.tsv.gz',
        'vogdb': 'a/fake/location/vogdb.hmm', 
        'pfam': 'tests/data/Pfam-A_subset.full.gz',
        'dbcan': None, 'peptidase': None, 'uniref': 'a/fake/location/uniref.fasta'}
    db_handler.filter_db_locs(low_mem_mode=True)
    assert db_handler.config.get('search_databases') == {
        'kofam_hmm': 'a/fake/location/kofam.hmm',
        'kofam_ko_list': 'a/fake/location/kofam_ko_list.tsv.gz',
        'pfam': 'tests/data/Pfam-A_subset.full.gz',
        'dbcan': None, 'peptidase': None}
    db_handler.config['search_databases'] = {
        'kegg': 'a/fake/location/kegg.fasta', 
        'kofam_hmm': 'a/fake/location/kofam.hmm',
        'kofam_ko_list': 'a/fake/location/kofam_ko_list.tsv.gz',
        'vogdb': 'a/fake/location/vogdb.hmm', 
        'pfam': 'tests/data/Pfam-A_subset.full.gz',
        'dbcan': None, 'peptidase': None}
    db_handler.filter_db_locs(use_vogdb=False)
    assert db_handler.config.get('search_databases') == {
        'kegg': 'a/fake/location/kegg.fasta', 
        'kofam_hmm': 'a/fake/location/kofam.hmm',
        'kofam_ko_list': 'a/fake/location/kofam_ko_list.tsv.gz',
        'pfam': 'tests/data/Pfam-A_subset.full.gz', 
        'dbcan': None, 
        'peptidase': None}
    db_handler.config['search_databases'] = {
        'kegg': 'a/fake/location/kegg.fasta', 
        'kofam_hmm': 'a/fake/location/kofam.hmm',
        'kofam_ko_list': 'a/fake/location/kofam_ko_list.tsv.gz',
        'vogdb': 'a/fake/location/vogdb.hmm', 
        'pfam': 'tests/data/Pfam-A_subset.full.gz',
        'dbcan': None, 'peptidase': None, 
        'uniref': 'a/fake/location/uniref.fasta'}
    db_handler.filter_db_locs(use_uniref=False)
    assert db_handler.config.get('search_databases') == {
        'kegg': 'a/fake/location/kegg.fasta', 
        'kofam_hmm': 'a/fake/location/kofam.hmm',
        'kofam_ko_list': 'a/fake/location/kofam_ko_list.tsv.gz',
        'vogdb': 'a/fake/location/vogdb.hmm', 
        'pfam': 'tests/data/Pfam-A_subset.full.gz',
        'dbcan': None, 'peptidase': None}
    # old version
    db_handler.config['search_databases'] = {
        'uniref': 'a fake locations', 
        'kegg': '/a/fake/loc', 'viral': ''}
    db_handler.filter_db_locs(master_list=('uniref', 'kegg'))
    assert db_handler.config.get('search_databases') == {
        'uniref': 'a fake locations', 
        'kegg': '/a/fake/loc'}
    # test
    db_handler.config['search_databases'] = {
        'uniref': 'a fake locations', 
        'kegg': '/a/fake/loc', 'viral': ''}
    with pytest.raises(ValueError):
        db_handler.filter_db_locs(low_mem_mode=True)
    # test
    db_handler.config['search_databases'] = {
        'kofam_hmm': '1', 'kofam_ko_list': '2', 
        'kegg': '/a/fake/loc', 'viral': ''}
    db_handler.filter_db_locs(low_mem_mode=True)
    assert db_handler.config.get('search_databases') == {
        'kofam_hmm': '1', 'kofam_ko_list': '2', 'viral': ''}
    # test
    db_handler.config['search_databases'] = {
        'uniref': 'a fake locations', 'kegg': '/a/fake/loc', 'viral': ''}
    db_handler.filter_db_locs(use_uniref=True)
    assert db_handler.config.get('search_databases') == {
        'uniref': 'a fake locations', 'kegg': '/a/fake/loc', 'viral': ''}
    # test
    db_handler.config['search_databases'] = {'uniref': 'a fake locations', 'kegg': '/a/fake/loc', 'viral': ''}
    db_handler.filter_db_locs(use_uniref=False)
    assert db_handler.config.get('search_databases') == {'kegg': '/a/fake/loc', 'viral': ''}


def test_make_header_dict_from_mmseqs_db():
    test_headers = DatabaseHandler.make_header_dict_from_mmseqs_db(path.join('tests', 'data', 'NC_001422.mmsdb'))
    assert test_headers == [{'id': 'NP_040705.1', 'description': 'NP_040705.1 B [Escherichia virus phiX174]'},
                            {'id': 'NP_040704.1', 'description': 'NP_040704.1 A* [Escherichia virus phiX174]'},
                            {'id': 'NP_040703.1', 'description': 'NP_040703.1 A [Escherichia virus phiX174]'},
                            {'id': 'NP_040713.1', 'description': 'NP_040713.1 H [Escherichia virus phiX174]'},
                            {'id': 'NP_040712.1', 'description': 'NP_040712.1 G [Escherichia virus phiX174]'},
                            {'id': 'NP_040711.1', 'description': 'NP_040711.1 F [Escherichia virus phiX174]'},
                            {'id': 'NP_040710.1', 'description': 'NP_040710.1 J [Escherichia virus phiX174]'},
                            {'id': 'NP_040709.1', 'description': 'NP_040709.1 E [Escherichia virus phiX174]'},
                            {'id': 'NP_040708.1', 'description': 'NP_040708.1 D [Escherichia virus phiX174]'},
                            {'id': 'NP_040707.1', 'description': 'NP_040707.1 C [Escherichia virus phiX174]'},
                            {'id': 'NP_040706.1', 'description': 'NP_040706.1 K [Escherichia virus phiX174]'}]


def test_process_pfam_descriptions():
    description_list = DatabaseHandler.process_pfam_descriptions(path.join('tests', 'data', 'Pfam-A_subset.hmm.dat.gz'))
    assert description_list == [{'id': 'PF10417.9', 'description': 'C-terminal domain of 1-Cys peroxiredoxin'},
                                {'id': 'PF12574.8', 'description': '120 KDa Rickettsia surface antigen'},
                                {'id': 'PF09847.9', 'description': 'Membrane protein of 12 TMs'},
                                {'id': 'PF00244.20', 'description': '14-3-3 protein'},
                                {'id': 'PF16998.5', 'description': '17 kDa outer membrane surface antigen'}]


def test_process_dbcan_descriptions():
    description_list = DatabaseHandler.process_dbcan_descriptions(
        path.join('tests', 'data', 'CAZyDB.07312019.fam-activities.subset.txt'),
        path.join('tests', 'data', 'example.fam.subfam.ec.txt'),
    )
    expected = pd.DataFrame([
        {'id': 'AA0', 'description': 'AA0', 'ec': ''},
        {'id': 'AA10', 
         'description': 'AA10 (formerly CBM33) proteins are copper-dependent '
                        'lytic polysaccharide monooxygenases (LPMOs); some '
                        'proteins have been shown to act on chitin, others on '
                        'cellulose; lytic cellulose monooxygenase '
                        '(C1-hydroxylating) (EC 1.14.99.54); lytic cellulose '
                        'monooxygenase (C4-dehydrogenating)(EC 1.14.99.56); '
                        'lytic chitin monooxygenase (EC 1.14.99.53)',
         'ec': '3.2.1.78,1.14.99.53,1.14.99.54'
         },
         {'id': 'AA11', 
          'description': 'AA11 proteins are copper-dependent lytic polysaccharide '
                         'monooxygenases (LPMOs); cleavage of chitin chains with '
                         'oxidation of C-1 has been demonstrated for a AA11 LPMO '
                         'from Aspergillus oryzae;',
          'ec': '1.14.99.55'
          },
         {'id': 'AA12', 
          'description': 'AA12 The pyrroloquinoline quinone-dependent oxidoreductase '
                         'activity was demonstrated for the CC1G_09525 protein of '
                         'Coprinopsis cinerea.',
          'ec': '1.14.99.56'}
          ])
    
    assert pd.DataFrame(description_list).equals(expected)


def test_process_vogdb_descriptions():
    description_list = DatabaseHandler.process_vogdb_descriptions(
        path.join('tests', 'data', 'vog_annotations_latest.subset.tsv.gz'))
    assert description_list == [{'id': 'VOG00001', 'description': 'sp|Q5UQJ2|YR863_MIMIV Putative ankyrin repeat '
                                                                  'protein R863; Xh'},
                                {'id': 'VOG00002', 'description': 'sp|Q9J4Z6|V244_FOWPN Putative ankyrin repeat '
                                                                  'protein FPV244; Xh'},
                                {'id': 'VOG00003', 'description': 'sp|Q91FD6|388R_IIV6 Putative MSV199 '
                                                                  'domain-containing protein 388R; Xu'},
                                {'id': 'VOG00004', 'description': 'sp|P51704|RPC1_BPHC1 Repressor protein CI; Xu'},
                                {'id': 'VOG00005', 'description': 'sp|P17766|POLG_PPVNA Genome polyprotein; XhXrXs'}]


def test_export_config(tmpdir):
    test_exported_config = path.join(tmpdir.mkdir('test_import_config'), 'test_CONFIG')
    export_config(test_exported_config)
    assert isinstance(json.loads(open(test_exported_config).read()), dict)
    assert isinstance(json.loads(open(test_exported_config).read())['search_databases'], dict)
    assert isinstance(json.loads(open(test_exported_config).read())['database_descriptions'], dict)
    assert isinstance(json.loads(open(test_exported_config).read())['dram_sheets'], dict)

def test_import_config(tmpdir):
    stash_conf = path.join(tmpdir.mkdir('stash'), 'stash')
    export_config(stash_conf)
    import_config(TEST_CONF_PATH)
    test_exported_config = path.join(tmpdir.mkdir('test_import_config'), 'test_CONFIG')
    export_config(test_exported_config)
    import_config(stash_conf)
    assert isinstance(json.loads(open(test_exported_config).read()), dict)
    assert isinstance(json.loads(open(test_exported_config).read())['search_databases'], dict)
    assert isinstance(json.loads(open(test_exported_config).read())['database_descriptions'], dict)
    assert isinstance(json.loads(open(test_exported_config).read())['dram_sheets'], dict)

def test_import_config_old(tmpdir):
    stash_conf = path.join(tmpdir.mkdir('stash'), 'stash')
    export_config(stash_conf)
    import_config(TEST_CONF_PATH_1_3)
    test_exported_config = path.join(tmpdir.mkdir('test_import_config_old'), 'test_CONFIG')
    export_config(test_exported_config)
    import_config(stash_conf)
    assert isinstance(json.loads(open(test_exported_config).read()), dict)
    assert isinstance(json.loads(open(test_exported_config).read())['search_databases'], dict)
    assert isinstance(json.loads(open(test_exported_config).read())['database_descriptions'], dict)
    assert isinstance(json.loads(open(test_exported_config).read())['dram_sheets'], dict)

def test_write_config(tmpdir, db_handler):
    stash_conf = path.join(tmpdir.mkdir('stash'), 'stash')
    export_config(stash_conf)
    test_exported_config = path.join(tmpdir.mkdir('test_write_config'), 'test_CONFIG')
    db_handler.write_config(test_exported_config)
    import_config(stash_conf)
    assert isinstance(json.loads(open(test_exported_config).read()), dict)
    assert isinstance(json.loads(open(test_exported_config).read())['search_databases'], dict)
    assert isinstance(json.loads(open(test_exported_config).read())['database_descriptions'], dict)
    assert isinstance(json.loads(open(test_exported_config).read())['dram_sheets'], dict)


@pytest.fixture()
def db_handler_1_3(logger):
    return DatabaseHandler(logger, TEST_CONF_PATH_1_3)


def test_construct_from_dram_pre_1_4_0(tmpdir, db_handler_1_3):
    test_exported_config = path.join(tmpdir.mkdir('test_export'), 'test_CONFIG')
    db_handler_1_3.write_config(test_exported_config)
    assert isinstance(json.loads(open(test_exported_config).read()), dict)
    assert isinstance(json.loads(open(test_exported_config).read())['search_databases'], dict)
    assert isinstance(json.loads(open(test_exported_config).read())['database_descriptions'], dict)
    assert isinstance(json.loads(open(test_exported_config).read())['dram_sheets'], dict)


def test_print_database_locations(capsys):
    print_database_locations()
    out, err = capsys.readouterr()
    assert 'Description db: ' in out
    assert err.endswith("Logging to console\n")

def test_enviroment_variable(logger):
    mypath:str = 'my_path_that-is_unique'
    environ['DRAM_CONFIG_LOCATION'] = mypath
    assert get_config_loc() == mypath
    del(environ['DRAM_CONFIG_LOCATION'])


def test_set_database_paths(logger):
    # first test that adding nothing doesn't change CONFIG
    db_handler = DatabaseHandler(logger)
    pretest_db_dict = db_handler.config.get('search_databases')
    db_handler.set_database_paths(write_config=False)
    assert pretest_db_dict == db_handler.config.get('search_databases')
    # test that adding something that doesn't exist throws error
    test_fake_database = path.join('tests', 'data', 'fake_database.mmsdb')
    with pytest.raises(ValueError):
        db_handler.set_database_paths(kegg_loc=test_fake_database, write_config=False)
    
    # test that adding something real is really added
    kegg_loc = path.join('tests', 'data', 'fake_gff.gff')
    db_handler.set_database_paths(kegg_loc=kegg_loc, write_config=False)
    assert db_handler.config.get('search_databases')['kegg'] == path.realpath(kegg_loc)


def test_clear_config(logger):
    db_handler = DatabaseHandler(logger)
    db_handler.clear_config(write_config=False)
    assert np.all([i is None for i in db_handler.config['search_databases'].values()])
    assert np.all([i is None for i in db_handler.config['database_descriptions'].values()])
    assert db_handler.config['description_db'] is None
