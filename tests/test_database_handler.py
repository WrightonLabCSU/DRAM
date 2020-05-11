import pytest
from os import path
from mag_annotator.database_setup import create_description_db, KeggDescription
from mag_annotator.database_handler import DatabaseHandler
from sqlalchemy.orm.exc import NoResultFound


@pytest.fixture()
def db_loc(tmpdir):
    db_loc = str(tmpdir.mkdir('test_db').join('test_db.sqlite'))
    create_description_db(db_loc)
    return db_loc


def test_create_description_db(db_loc):
    assert path.isfile(db_loc)


@pytest.fixture()
def db_handler(db_loc):
    return DatabaseHandler(db_loc)


def test_init_db_handler():
    with pytest.raises(ValueError):
        DatabaseHandler('my.fake.path.123.db.sqlite')


def test_add_descriptions_to_database(db_handler):
    kegg_entries = [{'id': 'K00001', 'description': 'The first KO'},
                    {'id': 'K00002', 'description': 'The second KO'}]
    db_handler.add_descriptions_to_database(kegg_entries, 'kegg_description')
    assert len(db_handler.session.query(KeggDescription).all()) == 2
    # now test that clear table works
    kegg_entry = [{'id': 'K00003', 'description': 'The third KO'}]
    db_handler.add_descriptions_to_database(kegg_entry, 'kegg_description', clear_table=True)
    assert len(db_handler.session.query(KeggDescription).all()) == 1


@pytest.fixture()
def db_w_entries(db_handler):
    db_handler.add_descriptions_to_database([{'id': 'K00001', 'description': 'The first KO'},
                                             {'id': 'K00002', 'description': 'The second KO'}], 'kegg_description',
                                            clear_table=True)
    return db_handler


def test_get_description(db_w_entries):
    description = db_w_entries.get_description('K00001', 'kegg_description')
    assert description == 'The first KO'
    with pytest.raises(NoResultFound):
        db_w_entries.get_description('K00003', 'kegg_description')


def test_get_descriptions(db_w_entries):
    description_dict = db_w_entries.get_descriptions(['K00001', 'K00002'], 'kegg_description')
    assert type(description_dict) is dict
    assert len(description_dict) == 2
    description_dict = db_w_entries.get_descriptions(['K00001', 'K00002', 'K00003'], 'kegg_description')
    assert type(description_dict) is dict
    assert len(description_dict) == 2
    description_dict = db_w_entries.get_descriptions(['K00003'], 'kegg_description')
    assert type(description_dict) is dict
    assert len(description_dict) == 0


def test_get_database_names(db_w_entries):
    names = db_w_entries.get_database_names()
    assert len(names) == 7
