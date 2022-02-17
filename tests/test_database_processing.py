import pytest
import os

from mag_annotator.database_processing import generate_modified_kegg_fasta, process_kegg_db, process_mmspro

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
    # assert os.path.isfile(os.path.join(processed_mmspro, 'fake.mmspro.idx'))
    assert os.path.isfile(os.path.join(processed_mmspro, 'fake.mmspro_h.index'))
