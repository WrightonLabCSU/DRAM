import pytest
import os
import logging

from mag_annotator.utils import setup_logger
from mag_annotator.database_processing import generate_modified_kegg_fasta, process_kegg, process_mmspro

@pytest.fixture()
def phix_proteins():
    return os.path.join('tests', 'data', 'NC_001422.faa')


@pytest.fixture()
def fake_gene_link_loc():
    return os.path.join('tests', 'data', 'fake_gene_ko.link')

@pytest.fixture()
def logger(tmpdir):
    logger = logging.getLogger('test_log')
    setup_logger(logger)
    return logger


def test_generate_modified_kegg_fasta(phix_proteins, fake_gene_link_loc):
    seqs = generate_modified_kegg_fasta(phix_proteins,
                                        fake_gene_link_loc)
    seqs_dict = {seq.metadata['id']: seq for seq in seqs}
    assert 'K14021' in seqs_dict['NP_040705.1'].metadata['description']
    assert 'K04464' in seqs_dict['NP_040708.1'].metadata['description']


def test_process_kegg(tmpdir, phix_proteins, fake_gene_link_loc, logger):
    processed_kegg_dir = tmpdir.mkdir('process_kegg_test')
    kegg = process_kegg(phix_proteins, processed_kegg_dir, logger, fake_gene_link_loc, download_date='today')
    assert os.path.isfile(kegg['kegg'])
    assert os.path.isfile('%s_h' % kegg['kegg'])


def test_process_mmspro(tmpdir, logger):
    processed_mmspro = tmpdir.mkdir('process_mmspro')
    process_mmspro(os.path.join('tests', 'data', 'Pfam-A_subset.full.gz'), processed_mmspro, logger, 'fake', 1, False)
    assert os.path.isfile(os.path.join(processed_mmspro, 'fake.mmspro'))
    assert os.path.isfile(os.path.join(processed_mmspro, 'fake.mmspro_h'))
    # assert os.path.isfile(os.path.join(processed_mmspro, 'fake.mmspro.idx'))
    assert os.path.isfile(os.path.join(processed_mmspro, 'fake.mmspro_h.index'))


# def test_down_and_process_camper(tmpdir, logger):
#     output_dir = tmpdir.mkdir('process_camper')
#     temporary = os.path.join(output_dir, 'temp')
#     os.mkdir(temporary)
#     loc= download_camper_tar_gz(temporary, logger)
#     process_camper(loc, output_dir, logger)

