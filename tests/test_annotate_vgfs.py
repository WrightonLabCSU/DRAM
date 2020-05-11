import pytest

import os
import pandas as pd
import numpy as np

from mag_annotator.annotate_vgfs import calculate_auxiliary_scores, get_virsorter_hits, get_overlap, get_gene_order, \
    is_transposon, get_metabolic_flags, get_amg_ids, get_virsorter_affi_contigs_name


@pytest.fixture()
def get_data_dir():
    return os.path.join(os.path.dirname(__file__), 'data')


def test_get_virsorter_hits(get_data_dir):
    virsorter_affi_contigs_tab = os.path.join(get_data_dir, 'VIRSorter_affi-contigs.tab')
    virsorter_hits = get_virsorter_hits(virsorter_affi_contigs_tab)
    assert virsorter_hits.shape == (14, 12)
    assert set(virsorter_hits['name']) == {'VIRSorter_DKWLOControl_V1KH_1_scaffold_4812-circular',
                                           'VIRSorter_DKWLPControl_V1KH_1_scaffold_5216'}
    assert virsorter_hits.loc['VIRSorter_DKWLPControl_V1KH_1_scaffold_5216-gene_3', 'viral_protein_cluster_category'] \
        == '1'


def test_get_overlap():
    f_overlap, r_overlap = get_overlap({'start_position': 0, 'end_position': 100},
                                       {'start_position': 101, 'end_position': 200})
    assert f_overlap == 0
    assert r_overlap == 0
    f_overlap, r_overlap = get_overlap({'start_position': 50, 'end_position': 150},
                                       {'start_position': 100, 'end_position': 200})
    assert f_overlap == .5
    assert r_overlap == .5


def test_is_transposon():
    assert not is_transposon(np.NaN)
    assert not is_transposon('ABCxyz')
    assert not is_transposon('PF11111.22')
    assert is_transposon('some things about it [PF01610.3]')
    assert is_transposon('some things about it [PF02371.3]; other things [PF00222.1]')


@pytest.fixture()
def gene_order1():
    return [('dram_gene1', 'virsorter_gene1', '-'),
            ('dram_gene2', 'virsorter_gene2', '4'),
            ('dram_gene3', 'virsorter_gene3', '-'),
            ('dram_gene4', 'virsorter_gene4', '1'),
            ('dram_gene5', 'virsorter_gene5', '-')]


@pytest.fixture()
def dram_genes():
    return pd.DataFrame([['dram_gene1', 0, 10],
                         ['dram_gene2', 12, 22],
                         ['dram_gene3', 24, 34],
                         ['dram_gene4', 36, 46],
                         ['dram_gene5', 48, 58]],
                        columns=['gene_name', 'start_position', 'end_position']).set_index('gene_name')


@pytest.fixture()
def virsorter_genes():
    return pd.DataFrame([['virsorter_gene1', 0, 10, '-'],
                         ['virsorter_gene2', 12, 22, '4'],
                         ['virsorter_gene3', 24, 34, '-'],
                         ['virsorter_gene4', 36, 46, '1'],
                         ['virsorter_gene5', 48, 58, '-']],
                        columns=['gene_name', 'start_position', 'end_position',
                                 'viral_protein_cluster_category']).set_index('gene_name')


def test_get_gene_order(dram_genes, virsorter_genes, gene_order1):
    test_gene_order1 = get_gene_order(dram_genes, virsorter_genes, min_overlap=.95)

    assert test_gene_order1 == gene_order1
    offset_gene_order = [('dram_gene1', None, None),
                         (None, 'virsorter_gene1', '-'),
                         ('dram_gene2', 'virsorter_gene2', '4'),
                         (None, 'virsorter_gene3', '-'),
                         ('dram_gene3', None, None),
                         ('dram_gene4', 'virsorter_gene4', '1'),
                         ('dram_gene5', 'virsorter_gene5', '-'),
                         (None, 'virsorter_gene6', '2'),
                         ('dram_gene6', None, None)]
    offset_dram_genes = pd.DataFrame([['dram_gene1', 0, 10],
                                      ['dram_gene2', 17, 27],
                                      ['dram_gene3', 30, 40],
                                      ['dram_gene4', 45, 55],
                                      ['dram_gene5', 56, 72],
                                      ['dram_gene6', 75, 300]],
                                     columns=['gene_name', 'start_position', 'end_position']).set_index('gene_name')
    offset_virsorter_genes = pd.DataFrame([['virsorter_gene1', 11, 15, '-'],
                                           ['virsorter_gene2', 17, 27, '4'],
                                           ['virsorter_gene3', 29, 34, '-'],
                                           ['virsorter_gene4', 45, 55, '1'],
                                           ['virsorter_gene5', 59, 69, '-'],
                                           ['virsorter_gene6', 75, 76, '2']],
                                          columns=['gene_name', 'start_position', 'end_position',
                                                   'viral_protein_cluster_category']).set_index('gene_name')
    test_gene_order2 = get_gene_order(offset_dram_genes, offset_virsorter_genes, min_overlap=.5)
    assert test_gene_order2 == offset_gene_order

    offset1_gene_order = [('dram_gene1', None, None),
                          (None, 'virsorter_gene1', '-'),
                          ('dram_gene2', 'virsorter_gene2', '4'),
                          (None, 'virsorter_gene3', '-'),
                          ('dram_gene3', None, None),
                          ('dram_gene4', 'virsorter_gene4', '1'),
                          ('dram_gene5', 'virsorter_gene5', '-'),
                          ('dram_gene6', None, None),
                          (None, 'virsorter_gene6', '2')]
    offset1_dram_genes = pd.DataFrame([['dram_gene1', 0, 10],
                                      ['dram_gene2', 17, 27],
                                      ['dram_gene3', 30, 40],
                                      ['dram_gene4', 45, 55],
                                      ['dram_gene5', 56, 72],
                                      ['dram_gene6', 75, 76]],
                                      columns=['gene_name', 'start_position', 'end_position']).set_index('gene_name')
    offset1_virsorter_genes = pd.DataFrame([['virsorter_gene1', 11, 15, '-'],
                                           ['virsorter_gene2', 17, 27, '4'],
                                           ['virsorter_gene3', 29, 34, '-'],
                                           ['virsorter_gene4', 45, 55, '1'],
                                           ['virsorter_gene5', 59, 69, '-'],
                                           ['virsorter_gene6', 75, 300, '2']],
                                           columns=['gene_name', 'start_position', 'end_position',
                                                    'viral_protein_cluster_category']).set_index('gene_name')
    test_gene_order3 = get_gene_order(offset1_dram_genes, offset1_virsorter_genes, min_overlap=.5)
    assert test_gene_order3 == offset1_gene_order


def test_get_auxilary_score(gene_order1):
    scores1 = calculate_auxiliary_scores(gene_order1)
    assert len(scores1) == 5
    assert scores1 == {'dram_gene1': 5,
                       'dram_gene2': 4,
                       'dram_gene3': 3,
                       'dram_gene4': 4,
                       'dram_gene5': 5}
    ones_and_twos_order = [('dram_gene1', 'virsorter_gene1', '-'),
                           ('dram_gene2', 'virsorter_gene2', '0'),
                           ('dram_gene3', 'virsorter_gene3', '4'),
                           ('dram_gene4', 'virsorter_gene4', '0'),
                           ('dram_gene5', 'virsorter_gene5', '4')]
    scores2 = calculate_auxiliary_scores(ones_and_twos_order)
    assert scores2 == {'dram_gene1': 5,
                       'dram_gene2': 4,
                       'dram_gene3': 1,
                       'dram_gene4': 2,
                       'dram_gene5': 5}

    five_no_others = [('dram_gene1', 'virsorter_gene1', '-'),
                      ('dram_gene2', 'virsorter_gene2', '-'),
                      ('dram_gene3', 'virsorter_gene3', '-'),
                      ('dram_gene4', 'virsorter_gene4', '-'),
                      ('dram_gene5', 'virsorter_gene5', '-')]
    scores3 = calculate_auxiliary_scores(five_no_others)
    assert scores3 == {'dram_gene1': 5,
                       'dram_gene2': 5,
                       'dram_gene3': 5,
                       'dram_gene4': 5,
                       'dram_gene5': 5}

    only_one = [('dram_gene1', 'virsorter_gene1', '-'),
                ('dram_gene2', 'virsorter_gene2', '-'),
                ('dram_gene3', 'virsorter_gene3', '0'),
                ('dram_gene4', 'virsorter_gene4', '-'),
                ('dram_gene5', 'virsorter_gene5', '-')]
    scores4 = calculate_auxiliary_scores(only_one)
    assert scores4 == {'dram_gene1': 5,
                       'dram_gene2': 4,
                       'dram_gene3': 4,
                       'dram_gene4': 4,
                       'dram_gene5': 5}


def test_get_metabolic_flags():
    annotations = pd.DataFrame([['scaffold_1', 101, 110, 'Xr', False, 'K00001'],
                                ['scaffold_1', 234, 423, None, False, 'K00002'],
                                ['scaffold_2', 52, 105, 'Xs', True, 'K00001'],
                                ['scaffold_3', 500, 582, None, False, 'K12345'],
                                ['scaffold_4', 101, 110, 'Xh', False, 'K12345'],
                                ['scaffold_4', 111, 120, 'Xh', False, 'K12345'],
                                ['scaffold_4', 121, 130, 'Xh', False, 'K12345'],
                                ['scaffold_4', 131, 139, 'Xh', False, 'K12345']],
                               index=['scaffold_1_1', 'scaffold_1_2', 'scaffold_2_1', 'scaffold_3_1', 'scaffold_4_1',
                                      'scaffold_4_2', 'scaffold_4_3', 'scaffold_4_4'],
                               columns=['scaffold', 'start_position', 'end_position', 'vogdb_categories',
                                        'is_transposon', 'kegg_id'])
    metabolic_genes = {'K12345'}
    amgs = {'K00001'}
    verified_amgs = {'K00001'}
    length_dict = {'scaffold_1': 250,
                   'scaffold_2': 150,
                   'scaffold_3': 600,
                   'scaffold_4': 1000}
    flags1 = get_metabolic_flags(annotations, metabolic_genes, amgs, verified_amgs, length_dict, 100)
    assert flags1 == {'scaffold_1_1': 'VMKE',
                      'scaffold_1_2': 'F',
                      'scaffold_2_1': 'VMKETF',
                      'scaffold_3_1': 'MF',
                      'scaffold_4_1': 'MB',
                      'scaffold_4_2': 'MB',
                      'scaffold_4_3': 'MB',
                      'scaffold_4_4': 'MB'}
                      # 'scaffold_4_2': 'MJB',
                      # 'scaffold_4_3': 'MJB',
                      # 'scaffold_4_4': 'MJB'}


def test_get_amg_ids():
    amg_frame = pd.DataFrame([['K00001', None, None],
                              [None, 'EC:2.2.1.2', 'PF0001']],
                             columns=['KO', 'EC', 'PFAM'])
    test_amg_ids = get_amg_ids(amg_frame)
    assert test_amg_ids == {'K00001', 'EC:2.2.1.2', 'PF0001'}


def test_get_virsorter_affi_contigs_name():
    assert get_virsorter_affi_contigs_name('scaffold_1-cat_6') == 'scaffold_1'
    assert get_virsorter_affi_contigs_name('scaffold_1_gene_1_gene_34-475-29373-cat_4') == 'scaffold_1'
    with pytest.raises(ValueError):
        get_virsorter_affi_contigs_name('my.fake.path.123.db.sqlite')
