import pytest

import os
import logging
import pandas as pd
import networkx as nx
import altair as alt

from mag_annotator.summarize_genomes import fill_genome_summary_frame, summarize_rrnas, \
    summarize_trnas, build_module_net, get_module_step_coverage, make_module_coverage_df, make_module_coverage_frame, \
    make_module_coverage_heatmap, pairwise, first_open_paren_is_all, split_into_steps, is_ko, make_module_network, \
    get_module_coverage, make_etc_coverage_df, make_etc_coverage_heatmap, make_functional_df, make_functional_heatmap, \
    fill_liquor_dfs, make_liquor_heatmap, make_liquor_df, make_genome_summary, write_summarized_genomes_to_xlsx,\
    get_phylum_and_most_specific, get_ids_from_annotations_all, get_ids_from_annotations_by_row
from mag_annotator.utils import get_ordered_uniques, setup_logger


def test_get_ordered_uniques():
    assert get_ordered_uniques([1, 2, 3]) == [1, 2, 3]
    assert get_ordered_uniques([1, 1, 2, 3]) == [1, 2, 3]
    assert get_ordered_uniques([1, 2, 1, 3]) == [1, 2, 3]


@pytest.fixture()
def logger(tmpdir):
    logger = logging.getLogger('test_log')
    setup_logger(logger)
    return logger


@pytest.fixture()
def annotations():
    return pd.DataFrame([['genome', 'K00001'],
                         ['genome', pd.np.NaN]],
                        index=['genome_scaffold_1_1', 'genome_scaffold_1_2'],
                        columns=['fasta', 'ko_id'])


@pytest.fixture()
def genome_summary_frame():
    return pd.DataFrame(pd.DataFrame([['K00001', 'description', 'module1', 'main', 'header1', 'subheader1'],
                                      ['K12345', 'description2', 'module2', 'main', 'header1', 'subheader1']],
                                     columns=['gene_id', 'gene_description', 'module', 'sheet', 'header', 'subheader']))


@pytest.fixture()
def summarized_genomes():
    return pd.DataFrame([['K00001', 'description', 'module1', 'main', 'header1', 'subheader1', 1],
                         ['K12345', 'description2', 'module2', 'main', 'header1', 'subheader1', 0]],
                        columns=['gene_id', 'gene_description', 'module', 'sheet', 'header', 'subheader', 'genome'])


def test_get_ids_from_row(logger):
    in_data = pd.concat([pd.DataFrame({'ko_id': 'K00001,K00003'}, index=['id_set1']),
                         pd.DataFrame({'kegg_hit': 'Some text and then [EC:0.0.0.0]; also [EC:1.1.1.1]'}, index=['id_set2']),
                         pd.DataFrame({'peptidase_family': 'ABC1;BCD2'}, index=['id_set3']),
                         pd.DataFrame({'cazy_best_hit': 'GH4'}, index=['id_set4'])
                        ])
    out_data = get_ids_from_annotations_by_row(in_data)
    assert out_data['id_set1'] == {'K00001', 'K00003'}
    assert out_data['id_set2'] == {'EC:0.0.0.0', 'EC:1.1.1.1'}
    assert out_data['id_set3'] == {'ABC1', 'BCD2'}
    assert out_data['id_set4'] == {'GH4'}


def test_fill_genome_summary_frame(annotations, genome_summary_frame, summarized_genomes, logger):
    test_frame = fill_genome_summary_frame(annotations, genome_summary_frame, 'fasta', logger)
    pd.testing.assert_frame_equal(test_frame, summarized_genomes)


def test_summarize_rrnas():
    fake_rrnas = pd.DataFrame([['fake_NC_001422.1', 990, 1000, 101.0, '+', '16S rRNA', pd.np.NaN]],
                              columns=['scaffold', 'begin', 'end', 'e-value', 'strand', 'type', 'note'])
    test_rrna_frame = summarize_rrnas(fake_rrnas, 'scaffold')
    print(test_rrna_frame)
    real_rrna_frame = pd.DataFrame([['5S rRNA', '5S ribosomal RNA gene', 'rRNA', 'rRNA', '', '', 0],
                                    ['16S rRNA', '16S ribosomal RNA gene', 'rRNA', 'rRNA', '', '', 1],
                                    ['23S rRNA', '23S ribosomal RNA gene', 'rRNA', 'rRNA', '', '', 0]],
                                   columns=['gene_id', 'gene_description', 'module', 'sheet', 'header', 'subheader',
                                            'fake_NC_001422.1'])
    pd.testing.assert_frame_equal(test_rrna_frame, real_rrna_frame)


@pytest.fixture()
def fake_processed_trnas():
    return pd.read_csv(os.path.join('tests', 'data', 'fake_trnas.tsv'), sep='\t')


def test_summarize_trnas(fake_processed_trnas):
    test_summarized_trnas = summarize_trnas(fake_processed_trnas)
    summarized_trnas = pd.DataFrame([['Gly (GCC)', 'Gly tRNA with GCC Codon', 'Gly tRNA', 'tRNA', 'tRNA', '', 1., 0.],
                                     ['Gly (TCC)', 'Gly tRNA with TCC Codon', 'Gly tRNA', 'tRNA', 'tRNA', '', 0., 1.],
                                     ['Leu (TAA)', 'Leu tRNA with TAA Codon', 'Leu tRNA', 'tRNA', 'tRNA', '', 1., 0.],
                                     ['Tyr, pseudo (GTA)', 'Tyr pseudo tRNA with GTA Codon', 'Tyr tRNA', 'tRNA', 'tRNA',
                                      '', 0, 1]],
                                    columns=['gene_id', 'gene_description', 'module', 'sheet', 'header', 'subheader',
                                             'Cytophaga_hutchinsonii_ATCC_33406', 'other_organismii'])
    pd.testing.assert_frame_equal(test_summarized_trnas, summarized_trnas)


def test_make_genome_summary(annotations, genome_summary_frame, summarized_genomes, logger):
    test_genome_summary = make_genome_summary(annotations, genome_summary_frame, logger)
    assert type(test_genome_summary) is pd.DataFrame
    pd.testing.assert_frame_equal(test_genome_summary, summarized_genomes)


def test_write_summarized_genomes_to_xlsx(summarized_genomes, tmpdir):
    output_file = os.path.join(tmpdir.mkdir('genome_summary_xlsx'), 'genome_summary.xlsx')
    write_summarized_genomes_to_xlsx(summarized_genomes, output_file)
    assert os.path.isfile(output_file)


@pytest.fixture()
def test_module_net():
    module_frame = pd.DataFrame([['0,0', 'M12345', 'a module name', 'K00001'],
                                 ['1,0', 'M12345', 'a module name', 'K00002'],
                                 ['2,0', 'M12345', 'a module name', 'K00003']],
                                columns=['path', 'module', 'module_name', 'ko'])
    test_module_net = build_module_net(module_frame)
    return test_module_net


def test_build_module_net(test_module_net):
    real_module_net = nx.DiGraph(num_steps=0, module_id='M12345', module_name='a module name')
    real_module_net.add_edges_from([('0,0', 'end_step_0'), ('end_step_0', '1,0'),
                                    ('1,0', 'end_step_1'), ('end_step_1', '2,0'),
                                    ('2,0', 'end_step_2')])
    assert nx.is_isomorphic(test_module_net, real_module_net)


def test_get_module_step_coverage(test_module_net):
    test_coverages1 = get_module_step_coverage(set([]), test_module_net)
    assert test_coverages1 == (3, 0, 0, [])
    test_coverages2 = get_module_step_coverage({'K00001', 'K00003'}, test_module_net)
    assert test_coverages2 == (3, 2, 2/3, ['K00001', 'K00003'])
    test_coverages2 = get_module_step_coverage({'K00001', 'K00003', 'K00002', 'K12345'}, test_module_net)
    assert test_coverages2 == (3, 3, 1, ['K00001', 'K00002', 'K00003'])


@pytest.fixture()
def test_annotations_df():
    return pd.DataFrame([['', 'scaffold_1'],
                         ['K12345', 'scaffold_1'],
                         ['K00001', 'scaffold_1']],
                        index=['gene_1', 'gene_2', 'gene_3'], columns=['ko_id', 'scaffold'])


def test_make_module_coverage_df(test_annotations_df, test_module_net):
    test_module_coverage_df = make_module_coverage_df(test_annotations_df, {'M12345': test_module_net})
    module_coverage_df = pd.DataFrame([['a module name', 3, 1, 1/3, 1, 'K00001', 'gene_3']],
                                      index=['M12345'],
                                      columns=['module_name', 'steps', 'steps_present', 'step_coverage', 'ko_count',
                                               'kos_present', 'genes_present'])
    pd.testing.assert_frame_equal(test_module_coverage_df, module_coverage_df)


@pytest.fixture()
def module_coverage_frame():
    return pd.DataFrame([['scaffold_1', 'M12345', 'a module name', 3, 1, 1 / 3, 1, 'K00001', 'gene_3']],
                        columns=['genome', 'module', 'module_name', 'steps', 'steps_present', 'step_coverage',
                                 'ko_count', 'kos_present', 'genes_present'])


def test_make_module_coverage_frame(test_annotations_df, test_module_net, module_coverage_frame):
    test_module_coverage_frame = make_module_coverage_frame(test_annotations_df, {'M12345': test_module_net},
                                                            groupby_column='scaffold')
    pd.testing.assert_frame_equal(test_module_coverage_frame, module_coverage_frame)


def test_make_module_coverage_heatmap(module_coverage_frame):
    hm = make_module_coverage_heatmap(module_coverage_frame)
    assert type(hm) is alt.Chart


def test_pairwise():
    assert list(pairwise([1, 2, 3])) == [(1, 2), (2, 3)]


def test_first_open_paren_is_all():
    assert first_open_paren_is_all('()')
    assert first_open_paren_is_all('(K1+K2-(K3+K5-K4))')
    assert not first_open_paren_is_all('(K1+K2-)K3+K5-K4)')


def test_split_into_steps():
    true_steps = ['K00330', 'K00331+K00332,K00331+K13378,K13380']
    assert split_into_steps('K00330+(K00331+K00332,K00331+K13378,K13380)', '+') == true_steps


def test_is_ko():
    assert is_ko('K00000')
    assert not is_ko('K1')


@pytest.fixture()
def module_network():
    network = nx.DiGraph()
    network.add_edges_from([('start', 'K00001'), ('K00001', 'K00002'), ('K00001', 'K00003'), ('K00002', 'K00004'),
                            ('K00003', 'K00013'), ('K00013', 'K00004')])
    return network


def test_make_module_network(module_network):
    test_network, _ = make_module_network('K00001+(K00002,K00003+K00013)+K00004')
    assert nx.is_isomorphic(test_network, module_network)


def test_get_module_coverage(module_network):
    module_network.add_edge('K00004', 'end')
    assert (3, 1, 1/3, {'K00001'}, {'K00002', 'K00004'}) == get_module_coverage(module_network, {'K00001'})
    assert (4, 1, 1/4, {'K00003'}, {'K00001', 'K00013', 'K00004'}) == get_module_coverage(module_network, {'K00003'})
    assert (3, 1, 1/3, {'K00002'}, {'K00001', 'K00004'}) == get_module_coverage(module_network, {'K00002', 'K99999'})


@pytest.fixture()
def etc_module_df():
    return pd.DataFrame([['K00001+(K00002,K00003+K00013)+K00004', 'Complex I', 'oxidoreductase', 'M00000']],
                        columns=['definition', 'complex', 'module_name', 'module_id'])


@pytest.fixture()
def etc_coverage_df():
    return pd.DataFrame([['M00000', 'oxidoreductase', 'I', 'scaffold_1', 3, 1, 1/3, 'K00001',
                          'K00002,K00004', 'Complex I: oxidoreductase']],
                        columns=['module_id', 'module_name', 'complex', 'genome', 'path_length',
                                 'path_length_coverage', 'percent_coverage', 'genes', 'missing_genes',
                                 'complex_module_name'])


def test_make_etc_coverage_df(test_annotations_df, etc_module_df, etc_coverage_df, logger):
    test_etc_coverage_df = make_etc_coverage_df(etc_module_df, test_annotations_df, 'scaffold')
    pd.testing.assert_frame_equal(test_etc_coverage_df, etc_coverage_df)


def test_make_etc_coverage_heatmap(etc_coverage_df):
    hm = make_etc_coverage_heatmap(etc_coverage_df)
    assert type(hm) is alt.HConcatChart


@pytest.fixture()
def function_heatmap_form():
    return pd.DataFrame([['Category1', 'SubCategory1', 'A function', 'K00001, K99999', 'A long function name', ''],
                         ['Category1', 'SubCategory1', 'A function', 'K00002', 'A second long function name', ''],
                         ['Category1', 'SubCategory2', 'B function', 'K12345', 'A longer function name', '']],
                        columns=['category', 'subcategory', 'function_name', 'function_ids', 'long_function_name',
                                 'gene_symbol'])


@pytest.fixture()
def functional_df():
    return pd.DataFrame([['Category1', 'SubCategory1', 'A function', 'K00001',
                          'A long function name; A second long function name', '', 'scaffold_1', False,
                          'Category1: A function'],
                         ['Category1', 'SubCategory2', 'B function', 'K12345', 'A longer function name', '',
                          'scaffold_1', True, 'Category1: B function']],
                        columns=['category', 'subcategory', 'function_name', 'function_ids', 'long_function_name',
                                 'gene_symbol', 'genome', 'present', 'category_function_name'])


def test_make_functional_df(test_annotations_df, function_heatmap_form, functional_df):
    test_functional_df = make_functional_df(test_annotations_df, function_heatmap_form, logger, 'scaffold')
    pd.testing.assert_frame_equal(test_functional_df, functional_df)


def test_make_functional_heatmap(functional_df):
    hm = make_functional_heatmap(functional_df)
    assert type(hm) is alt.HConcatChart


# TODO: actually test that the frames are correct, already done above
def test_fill_liquor_dfs(test_annotations_df, test_module_net, etc_module_df, function_heatmap_form):
    module_nets = {'M12345': test_module_net}
    liquor_dfs = fill_liquor_dfs(test_annotations_df, module_nets, etc_module_df, function_heatmap_form, logger, 'scaffold')
    assert len(liquor_dfs) == 3
    assert type(liquor_dfs[0]) is pd.DataFrame
    assert type(liquor_dfs[1]) is pd.DataFrame
    assert type(liquor_dfs[2]) is pd.DataFrame


def test_make_liquor_heatmap(module_coverage_frame, etc_coverage_df, functional_df):
    hm = make_liquor_heatmap(module_coverage_frame, etc_coverage_df, functional_df)
    assert type(hm) is alt.HConcatChart


def test_make_liquor_df(module_coverage_frame, etc_coverage_df, functional_df):
    liquor_df = pd.DataFrame([[1/3, 1/3, False, True]],
                             index=pd.Index(['scaffold_1'], name='genome'),
                             columns=['a module name', 'Complex I: oxidoreductase', 'Category1: A function',
                                      'Category1: B function'])
    test_liquor_df = make_liquor_df(module_coverage_frame, etc_coverage_df, functional_df)
    pd.testing.assert_frame_equal(test_liquor_df, liquor_df)


def test_get_phylum_and_most_specific():
    assert get_phylum_and_most_specific('d__Bacteria;p__Bacteroidota;c__;o__;f__;g__;s__') == \
           'p__Bacteroidota;c__'
    assert get_phylum_and_most_specific('d__Archaea;p__;c__;o__;f__;g__;s__') == 'd__Archaea;p__'
    assert get_phylum_and_most_specific('d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;'
                                        'g__Alistipes;s__') == 'p__Bacteroidota;g__Alistipes'
    assert get_phylum_and_most_specific('d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;'
                                        'g__Enterococcus_D;s__Enterococcus_D gallinarum') == \
        'p__Firmicutes;s__Enterococcus_D gallinarum'
