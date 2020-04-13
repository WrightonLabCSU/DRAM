import pytest

import os
import json

from mag_annotator.utils import get_database_locs
from mag_annotator.database_processing import check_file_exists, make_header_dict_from_mmseqs_db, \
    generate_modified_kegg_fasta, process_kegg_db, process_mmspro, process_pfam_descriptions, \
    process_dbcan_descriptions, process_vogdb_descriptions, check_exists_and_add_to_location_dict, export_config, \
    print_database_locations, set_database_paths


def test_check_file_exists():
    assert check_file_exists(None)
    assert check_file_exists(os.path.join('tests', 'data', 'fake_gff.gff'))
    with pytest.raises(ValueError):
        check_file_exists(os.path.join('tests', 'data', 'a_fake_nonexistent_file.txt'))


def test_make_header_dict_from_mmseqs_db():
    test_headers = make_header_dict_from_mmseqs_db(os.path.join('tests', 'data', 'NC_001422.mmsdb'))
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


@pytest.fixture()
def phix_proteins():
    return os.path.join('tests', 'data', 'NC_001422.faa')


@pytest.fixture()
def fake_gene_link_loc():
    return os.path.join('tests', 'data', 'fake_gene_ko.link')


def test_generate_modified_kegg_fasta(phix_proteins, fake_gene_link_loc):
    seqs = generate_modified_kegg_fasta(phix_proteins,
                                        fake_gene_link_loc)
    seqs_dict = {seq.metadata['id']: seq for seq in seqs}
    assert 'K14021' in seqs_dict['NP_040705.1'].metadata['description']
    assert 'K04464' in seqs_dict['NP_040708.1'].metadata['description']


def test_process_kegg_db(tmpdir, phix_proteins, fake_gene_link_loc):
    processed_kegg_dir = tmpdir.mkdir('process_kegg_test')
    kegg_db = process_kegg_db(processed_kegg_dir, phix_proteins, fake_gene_link_loc, download_date='today')
    assert os.path.isfile(kegg_db)
    assert os.path.isfile('%s_h' % kegg_db)


def test_process_mmspro(tmpdir):
    processed_mmspro = tmpdir.mkdir('process_mmspro')
    process_mmspro(os.path.join('tests', 'data', 'Pfam-A_subset.full.gz'), processed_mmspro, 'fake', 1, False)
    assert os.path.isfile(os.path.join(processed_mmspro, 'fake.mmspro'))
    assert os.path.isfile(os.path.join(processed_mmspro, 'fake.mmspro_h'))
    assert os.path.isfile(os.path.join(processed_mmspro, 'fake.mmspro.idx'))
    assert os.path.isfile(os.path.join(processed_mmspro, 'fake.mmspro_h.index'))


def test_process_pfam_descriptions():
    description_list = process_pfam_descriptions(os.path.join('tests', 'data', 'Pfam-A_subset.hmm.dat.gz'))
    assert description_list == [{'id': 'PF10417.9', 'description': 'C-terminal domain of 1-Cys peroxiredoxin'},
                                {'id': 'PF12574.8', 'description': '120 KDa Rickettsia surface antigen'},
                                {'id': 'PF09847.9', 'description': 'Membrane protein of 12 TMs'},
                                {'id': 'PF00244.20', 'description': '14-3-3 protein'},
                                {'id': 'PF16998.5', 'description': '17 kDa outer membrane surface antigen'}]


def test_process_dbcan_descriptions():
    description_list = process_dbcan_descriptions(os.path.join('tests', 'data', 'CAZyDB.07312019.fam-activities.subset.txt'))
    assert description_list == [{'id': 'AA0', 'description': 'AA0'},
                                {'id': 'AA10', 'description': 'AA10 (formerly CBM33) proteins are copper-dependent '
                                                              'lytic polysaccharide monooxygenases (LPMOs); some '
                                                              'proteins have been shown to act on chitin, others on '
                                                              'cellulose; lytic cellulose monooxygenase '
                                                              '(C1-hydroxylating) (EC 1.14.99.54); lytic cellulose '
                                                              'monooxygenase (C4-dehydrogenating)(EC 1.14.99.56); '
                                                              'lytic chitin monooxygenase (EC 1.14.99.53)'},
                                {'id': 'AA11', 'description': 'AA11 proteins are copper-dependent lytic polysaccharide '
                                                              'monooxygenases (LPMOs); cleavage of chitin chains with '
                                                              'oxidation of C-1 has been demonstrated for a AA11 LPMO '
                                                              'from Aspergillus oryzae;'},
                                {'id': 'AA12', 'description': 'AA12 The pyrroloquinoline quinone-dependent oxidoreductase '
                                                              'activity was demonstrated for the CC1G_09525 protein of '
                                                              'Coprinopsis cinerea.'}]


def test_process_vogdb_descriptions():
    description_list = process_vogdb_descriptions(os.path.join('tests', 'data', 'vog_annotations_latest.subset.tsv.gz'))
    assert description_list == [{'id': 'VOG00001', 'description': 'sp|Q5UQJ2|YR863_MIMIV Putative ankyrin repeat '
                                                                  'protein R863; Xh'},
                                {'id': 'VOG00002', 'description': 'sp|Q9J4Z6|V244_FOWPN Putative ankyrin repeat '
                                                                  'protein FPV244; Xh'},
                                {'id': 'VOG00003', 'description': 'sp|Q91FD6|388R_IIV6 Putative MSV199 '
                                                                  'domain-containing protein 388R; Xu'},
                                {'id': 'VOG00004', 'description': 'sp|P51704|RPC1_BPHC1 Repressor protein CI; Xu'},
                                {'id': 'VOG00005', 'description': 'sp|P17766|POLG_PPVNA Genome polyprotein; XhXrXs'}]


def test_check_exists_and_add_to_location_dict(phix_proteins):
    test_dict1 = check_exists_and_add_to_location_dict(phix_proteins, 'fake', {})
    assert test_dict1 == {'fake': os.path.realpath(phix_proteins)}
    test_dict2 = check_exists_and_add_to_location_dict(None, 'fake', {})
    assert test_dict2 == {'fake': None}
    test_dict3 = check_exists_and_add_to_location_dict(None, 'fake', {'fake': 'abc123'})
    assert test_dict3 == {'fake': 'abc123'}


def test_export_config(capsys, tmpdir):
    export_config()
    out, err = capsys.readouterr()
    assert len(err) == 0
    assert 'description_db' in out
    assert out.startswith('{')
    assert out.strip().endswith('}')

    test_exported_config = os.path.join(tmpdir.mkdir('test_export'), 'test_CONFIG')
    export_config(test_exported_config)
    assert type(json.loads(open(test_exported_config).read())) is dict


def test_print_database_locations(capsys):
    print_database_locations()
    out, err = capsys.readouterr()
    assert 'Description db: ' in out
    assert len(err) == 0


def test_set_database_paths(tmpdir):
    test_config_dir = tmpdir.mkdir('test_config')
    # first test that adding nothing doesn't change CONFIG
    test_config = os.path.join(test_config_dir, 'CONFIG')
    pretest_db_dict = get_database_locs()
    set_database_paths(config_loc=test_config)
    test_db_dict = get_database_locs(test_config)
    assert type(test_db_dict) is dict
    assert pretest_db_dict == test_db_dict
    # test that adding something that doesn't exist throws error
    test_fake_database = os.path.join(test_config_dir, 'fake_database.mmsdb')
    with pytest.raises(ValueError):
        set_database_paths(kegg_db_loc=test_fake_database)
    # test that adding something real is really added
    kegg_loc = os.path.join('tests', 'data', 'fake_gff.gff')
    set_database_paths(kegg_db_loc=kegg_loc, config_loc=test_config)
    test_db_dict = get_database_locs(test_config)
    assert test_db_dict['kegg'] == os.path.realpath(kegg_loc)
    # test that adding something with use_current_locs False works
    set_database_paths(kegg_db_loc=kegg_loc, config_loc=test_config, use_current_locs=False)
    test_db_dict = get_database_locs(test_config)
    assert test_db_dict['kegg'] == os.path.realpath(kegg_loc)
    assert test_db_dict['description_db'] is None
