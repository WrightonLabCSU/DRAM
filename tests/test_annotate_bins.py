import pytest

import os

import pandas as pd

from mag_annotator.utils import make_mmseqs_db
from mag_annotator.annotate_bins import filter_fasta, run_prodigal, get_best_hits,\
    get_reciprocal_best_hits, process_reciprocal_best_hits, get_kegg_description, get_uniref_description, get_sig,\
    get_unannotated, assign_grades


@pytest.fixture()
def fasta_loc():
    return os.path.join('tests', 'data', 'NC_001422.fasta')


def test_filter_fasta(fasta_loc, tmpdir):
    filtered_seq_5000 = filter_fasta(fasta_loc, min_len=5000)
    assert len(filtered_seq_5000) == 1
    assert len(filtered_seq_5000[0]) == 5386
    filtered_seq_6000 = filter_fasta(fasta_loc, min_len=6000)
    assert len(filtered_seq_6000) == 0
    filt_fasta = tmpdir.mkdir('test_filt').join('filtered_fasta.fasta')
    filter_fasta(fasta_loc, min_len=5000, output_loc=str(filt_fasta))
    assert os.path.isfile(filt_fasta)


@pytest.fixture()
def prodigal_dir(fasta_loc, tmpdir):
    prodigal_output = tmpdir.mkdir('prodigal_output')
    gff, fna, faa = run_prodigal(fasta_loc, str(prodigal_output))
    return gff, fna, faa


@pytest.fixture()
def prodigal_gff(prodigal_dir):
    return prodigal_dir[0]


@pytest.fixture()
def prodigal_fna(prodigal_dir):
    return prodigal_dir[1]


@pytest.fixture()
def prodigal_faa(prodigal_dir):
    return prodigal_dir[2]


def test_run_prodigal(prodigal_gff, prodigal_fna, prodigal_faa):
    assert os.path.isfile(prodigal_gff)
    assert os.path.isfile(prodigal_fna)
    assert os.path.isfile(prodigal_faa)


@pytest.fixture()
def mmseqs_db_dir(tmpdir):
    output_loc = tmpdir.mkdir('make_mmseqs_db_test')
    return output_loc


@pytest.fixture()
def mmseqs_db(prodigal_faa, mmseqs_db_dir):
    output_file = str(mmseqs_db_dir.join('mmseqs_db.mmsdb'))
    make_mmseqs_db(prodigal_faa, output_file, True, 1)
    return output_file


def test_make_mmseqs_db(mmseqs_db):
    assert os.path.isfile(mmseqs_db)


@pytest.fixture()
def phix_proteins():
    return os.path.join('tests', 'data', 'NC_001422.faa')


@pytest.fixture()
def target_mmseqs_db(mmseqs_db_dir, phix_proteins):
    output_file = str(mmseqs_db_dir.join('target.mmsdb'))
    make_mmseqs_db(phix_proteins, output_file, True, 1)
    return output_file


@pytest.fixture()
def best_hits_loc(mmseqs_db, target_mmseqs_db, mmseqs_db_dir):
    best_hits_loc = get_best_hits(mmseqs_db, target_mmseqs_db, mmseqs_db_dir, threads=1, verbose=False)
    return best_hits_loc


def test_get_best_hits(best_hits_loc):
    assert os.path.isfile(best_hits_loc)


@pytest.fixture()
def reverse_best_hits_loc(best_hits_loc, mmseqs_db, target_mmseqs_db, mmseqs_db_dir):
    reverse_best_hits_loc = get_reciprocal_best_hits(mmseqs_db, target_mmseqs_db, mmseqs_db_dir, threads=1,
                                                     verbose=False)
    return reverse_best_hits_loc


def test_get_reciprocal_best_hits(reverse_best_hits_loc):
    assert os.path.isfile(reverse_best_hits_loc)


@pytest.fixture()
def processed_hits():
    forward = os.path.join('tests', 'data', 'query_target_hits.b6')
    reverse = os.path.join('tests', 'data', 'target_query_hits.b6')
    processed_hits = process_reciprocal_best_hits(forward, reverse)
    return processed_hits


def test_process_reciprocal_best_hits(processed_hits):
    assert processed_hits.shape == (7, 5)
    assert set(processed_hits.loc[processed_hits.target_RBH].index) == {'NC_001422.1_5', 'NC_001422.1_4',
                                                                        'NC_001422.1_7', 'NC_001422.1_6'}


def test_get_kegg_description():
    header_dict = {'aad:TC41_2367': 'aad:TC41_2367  ABC-type molybdate transport system periplasmic component-like protein',
                   'aar:Acear_0854': 'aar:Acear_0854  hypothetical protein; K05810 conserved hypothetical protein',
                   'aar:Acear_1520': 'aar:Acear_1520  hypothetical protein'}
    kegg_hits_data = [['aad:TC41_2367', 10e-5],
                      ['aar:Acear_0854', 10e-6],
                      ['aar:Acear_1520', 10e-10]]
    kegg_hits = pd.DataFrame(kegg_hits_data, index=['gene1', 'gene2', 'gene3'], columns=['kegg_hit', 'eVal'])
    kegg_hits_add_description = get_kegg_description(kegg_hits, header_dict)
    assert kegg_hits_add_description.shape == (3, 3)
    assert kegg_hits_add_description.loc['gene2', 'kegg_id'] == 'K05810'
    assert kegg_hits_add_description.loc['gene1', 'kegg_id'] == ''


def test_get_uniref_description():
    header_dict = {'UniRef90_A0A139CGD2': 'UniRef90_A0A139CGD2 Phosphate transport system permease protein PstA n=2 '
                                          'Tax=Candidatus Frackibacter TaxID=2017975 RepID=A0A139CGD2_9FIRM',
                   'UniRef90_2642661139': 'UniRef90_2642661139 Ga0073286_10147 conserved hypothetical protein n=1 '
                                          'Tax=Fuchsiella alkaliacetigena WG11 RepID=Ga0073286_10147'}
    uniref_hits_data = [['UniRef90_A0A139CGD2', 10e-5],
                        ['UniRef90_2642661139', 10e-20]]
    uniref_hits = pd.DataFrame(uniref_hits_data, index=['gene1', 'gene2'], columns=['uniref_hit', 'eVal'])
    uniref_hits_add_description = get_uniref_description(uniref_hits, header_dict)
    assert uniref_hits_add_description.shape == (2, 4)
    assert uniref_hits_add_description.loc['gene1', 'uniref_id'] == 'A0A139CGD2_9FIRM'
    assert uniref_hits_add_description.loc['gene1', 'uniref_taxonomy'] == 'Candidatus Frackibacter'


def test_get_viral_description():
    header_dict = {'YP_009015653.1': 'YP_009015653.1 gp350 [Bacillus virus G]',
                   'NP_077550.1': 'NP_077550.1 EsV-1-65 [Ectocarpus siliculosus virus 1]'}
    viral_hits_data = [['NP_077550.1'],
                       ['YP_009015653.1']]
    viral_hits = pd.DataFrame(viral_hits_data, index=['gene1', 'gene2'], columns=['viral_hit'])
    viral_hits_add_description = get_viral_description(viral_hits, header_dict)
    assert viral_hits_add_description.shape == (2, 2)
    assert viral_hits_add_description.loc['gene1', 'viral_id'] == 'NP_077550.1'


def test_run_mmseqs_pfam():
    pass


def test_get_sig():
    assert not get_sig(87, 1)
    assert get_sig(87, 1e-10)
    assert not get_sig(79, 1e-3)
    assert get_sig(79, 1e-19)


def test_run_hmmscan_dbcan():
    pass


def test_get_scaffold_and_gene():
    scaffold_gene_names = ['scaffold_0_1',
                           'scaffold_0_12',
                           'scaffold_3_13',
                           'scaffold_0_14',
                           'scaffold_99_15',
                           'scaffold_10_16']

    scaffold_gene_df = get_scaffold_and_gene(scaffold_gene_names)
    assert scaffold_gene_df.shape == (6, 2)
    assert scaffold_gene_df.loc['scaffold_0_12', 'scaffold'] == 'scaffold_0'
    assert scaffold_gene_df.loc['scaffold_0_12', 'gene_position'] == '12'
    assert scaffold_gene_df.loc['scaffold_99_15', 'scaffold'] == 'scaffold_99'
    assert scaffold_gene_df.loc['scaffold_99_15', 'gene_position'] == '15'


def test_get_unannotated(phix_proteins):
    annotated_genes = ['NP_040704.1', 'NP_040703.1', 'NP_040713.1', 'NP_040712.1', 'NP_040711.1', 'NP_040710.1',
                       'NP_040709.1', 'NP_040707.1']
    unannotated_genes = ['NP_040705.1', 'NP_040706.1', 'NP_040708.1']
    test_unannotated_genes = get_unannotated(phix_proteins, annotated_genes)
    assert set(unannotated_genes) == set(test_unannotated_genes)


def test_assign_grades():
    annotations_data = [[True, 'K00001', False, 'KOER09234OK', ['PF00001']],
                        [False, 'K00002', True, 'KLODKJFSO234KL', ['PF01234']],
                        [False, 'K00003', False, 'EIORWU234KLKDS', pd.np.NaN],
                        [False, pd.np.NaN, False, pd.np.NaN, pd.np.NaN]]
    annotations = pd.DataFrame(annotations_data, index=['gene1', 'gene2', 'gene3', 'gene4'],
                               columns=['kegg_RBH', 'kegg_hit', 'uniref_RBH', 'uniref_hit', 'pfam_hits'])
    test_grades = assign_grades(annotations)
    assert test_grades.loc['gene1'] == 'A'
    assert test_grades.loc['gene2'] == 'B'
    assert test_grades.loc['gene3'] == 'C'
    assert test_grades.loc['gene4'] == 'E'
