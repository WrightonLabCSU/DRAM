import pytest

from os import path

from mag_annotator.database_setup import create_description_db


@pytest.fixture()
def db_loc(tmpdir):
    db_loc = str(tmpdir.mkdir('test_db').join('test_db.sqlite'))
    create_description_db(db_loc)
    return db_loc


def test_create_description_db(db_loc):
    assert path.isfile(db_loc)
