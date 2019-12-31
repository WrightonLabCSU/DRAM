import pytest
from mag_annotator.annotate_vgfs import calculate_auxiliary_scores

@pytest.fixture()
def gene_order1():
    return (('dram_gene1', 'virsorter_gene1', '-'),
            ('dram_gene2', 'virsorter_gene2', '4'),
            ('dram_gene3', 'virsorter_gene3', '-'),
            ('dram_gene4', 'virsorter_gene4', '1'),
            ('dram_gene5', 'virsorter_gene5', '-'))


def test_get_auxilary_score(gene_order1):
    scores1 = calculate_auxiliary_scores(gene_order1)
    assert len(scores1) == 5
    assert scores1 == {'dram_gene1': 5,
                       'dram_gene2': 4,
                       'dram_gene3': 3,
                       'dram_gene4': 4,
                       'dram_gene5': 5}
