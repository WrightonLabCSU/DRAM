import pytest

import os
import pandas as pd

from mag_annotator.summarize_genomes import get_ordered_uniques, fill_genome_summary_frame, summarize_rrnas


def test_get_ordered_uniques():
    assert get_ordered_uniques([1, 2, 3]) == [1, 2, 3]
    assert get_ordered_uniques([1, 1, 2, 3]) == [1, 2, 3]
    assert get_ordered_uniques([1, 2, 1, 3]) == [1, 2, 3]


@pytest.fixture()
def annotations():
    return pd.DataFrame([['genome', 'K00001'],
                         ['genome', pd.np.NaN]],
                        index=['genome_scaffold_1_1', 'genome_scaffold_1_2'],
                        columns=['fasta', 'kegg_id'])


@pytest.fixture()
def genome_summary_frame():
    return pd.DataFrame([['K00001'],
                         ['K12345']],
                        columns=['gene_id'])


def test_fill_genome_summary_frame(annotations, genome_summary_frame):
    test_frame = fill_genome_summary_frame(annotations, genome_summary_frame, 'fasta')
    pd.testing.assert_frame_equal(test_frame, pd.DataFrame([['K00001', 1],
                                                            ['K12345', 0]],
                                                            columns=['gene_id', 'genome']))


def test_summary_rrnas():
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
