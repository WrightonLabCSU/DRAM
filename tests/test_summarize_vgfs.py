import pytest

import pandas as pd
import altair as alt

from mag_annotator.summarize_vgfs import filter_to_amgs, get_strand_switches, get_ids_from_row, make_viral_distillate, \
    make_vgf_order, make_amg_count_column, make_viral_functional_df, make_viral_functional_heatmap


@pytest.fixture()
def annotations():
    return pd.DataFrame([['scaffold_1', 'K99999', 'VTF', 1],
                         ['scaffold_1', 'K12345', 'M', 4],
                         ['scaffold_1', 'K00001', 'MFTJ', 3],
                         ['scaffold_1', 'K11111', 'MKE', 2]],
                        index=['gene_1', 'gene_2', 'gene_3', 'gene_4'],
                        columns=['scaffold', 'kegg_id', 'amg_flags', 'auxiliary_score'])


@pytest.fixture()
def pamgs():
    return pd.DataFrame([['scaffold_1', 'K00001', 'MFTJ', 3],
                         ['scaffold_1', 'K11111', 'MKE', 2]],
                        index=['gene_3', 'gene_4'],
                        columns=['scaffold', 'kegg_id', 'amg_flags', 'auxiliary_score'])


def test_filter_to_amgs(annotations, pamgs):
    test_pamgs = filter_to_amgs(annotations, remove_transposons=False, max_aux=3)
    pd.testing.assert_frame_equal(test_pamgs, pamgs)


def test_get_strand_switches():
    assert get_strand_switches(['+', '-', '+']) == 2
    assert get_strand_switches(['+', '+', '+']) == 0
    assert get_strand_switches(['-']) == 0
    assert get_strand_switches(['+', '+', '-', '-']) == 1


def test_get_ids_from_row():
    id_set1 = get_ids_from_row(pd.Series({'kegg_id': 'K00001,K00003'}))
    assert id_set1 == {'K00001', 'K00003'}
    id_set2 = get_ids_from_row(pd.Series({'kegg_hit': 'Some text and then [EC:0.0.0.0]; also [EC:1.1.1.1]'}))
    assert id_set2 == {'EC:0.0.0.0', 'EC:1.1.1.1'}
    id_set3 = get_ids_from_row(pd.Series({'peptidase_family': 'ABC1;BCD2'}))
    assert id_set3 == {'ABC1', 'BCD2'}
    id_set4 = get_ids_from_row(pd.Series({'cazy_hits': 'GH4 some things;GT6 other things'}))
    assert id_set4 == {'GH4', 'GT6'}


@pytest.fixture()
def genome_summary_frame():
    return pd.DataFrame(pd.DataFrame([['description', 'module1', 'main', 'header1', 'subheader1'],
                                      ['description2', 'module2', 'main', 'header1', 'subheader1'],
                                      ['description3', 'module3', 'other', 'header1', 'subheader1']],
                                     index=['K00001', 'K12345', 'K00001'],
                                     columns=['gene_description', 'module', 'sheet', 'header', 'subheader']))


def test_make_viral_distillate(pamgs, genome_summary_frame):
    test_viral_distillate = make_viral_distillate(pamgs, genome_summary_frame)
    viral_distillate = pd.DataFrame([['gene_3', 'scaffold_1', 'K00001', 'description', 'main', 'header1', 'subheader1',
                                      'module1', 3, 'MFTJ'],
                                     ['gene_3', 'scaffold_1', 'K00001', 'description3', 'other', 'header1',
                                      'subheader1', 'module3', 3, 'MFTJ'],
                                     ['gene_4', 'scaffold_1', '', '', '', '', '', '', 2, 'MKE']],
                                    columns=['gene', 'scaffold', 'gene_id', 'gene_description', 'category', 'header',
                                             'subheader', 'module', 'auxiliary_score', 'amg_flags'])
    pd.testing.assert_frame_equal(test_viral_distillate, viral_distillate)


def test_make_vgf_order(pamgs):
    assert make_vgf_order(pamgs) == ['scaffold_1']


def test_make_amg_count_column(pamgs):
    assert type(make_amg_count_column(pamgs)) is alt.Chart


@pytest.fixture()
def viral_functional_df():
    return pd.DataFrame([['main', 'module1', 'gene_3', 'K00001', 'scaffold_1', True],
                         ['other', 'module3', 'gene_3', 'K00001', 'scaffold_1', True]],
                        columns=['Category', 'Function', 'AMG Genes', 'Genes Present', 'Contig Name',
                                 'Present in Contig'])


def test_make_viral_functional_df(pamgs, genome_summary_frame, viral_functional_df):
    test_viral_function_df = make_viral_functional_df(pamgs, genome_summary_frame, 'scaffold')
    pd.testing.assert_frame_equal(test_viral_function_df, viral_functional_df)


def test_make_viral_functional_heatmap(viral_functional_df):
    assert type(make_viral_functional_heatmap(viral_functional_df)) is alt.HConcatChart
