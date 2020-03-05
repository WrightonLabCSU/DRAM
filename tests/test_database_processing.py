import pytest

import os

from mag_annotator.database_processing import check_file_exists, make_header_dict_from_mmseqs_db, remove_prefix


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


def test_remove_prefix():
    assert remove_prefix('prefix', 'pre') == 'fix'
    assert remove_prefix('postfix', 'pre') == 'postfix'


def test_generate_modified_kegg_fasta():
    pass
