from os import path
from warnings import warn

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from mag_annotator.database_setup import TABLE_NAME_TO_CLASS_DICT

# TODO: store all sequence db locs within database handler class


class DatabaseHandler:
    def __init__(self, database_loc):
        if not path.exists(database_loc):
            raise ValueError('Database does not exist at this path')
        engine = create_engine('sqlite:///%s' % database_loc)
        DBSession = sessionmaker(bind=engine)
        self.session = DBSession()

    # functions for adding descriptions to tables
    def add_descriptions_to_database(self, description_dict, db_name, clear_table=True):
        description_class = TABLE_NAME_TO_CLASS_DICT[db_name]
        if clear_table:
            self.session.query(description_class).delete()
        for annotation_id, description in description_dict.items():  # TODO: Make this a bulk operation
            self.session.add(description_class(id=annotation_id, description=description))
        self.session.commit()

    # functions for getting descriptions from tables
    def get_description(self, annotation_id, db_name):
        return self.session.query(TABLE_NAME_TO_CLASS_DICT[db_name]).filter_by(id=annotation_id).one().description

    def get_descriptions(self, ids, db_name):
        description_class = TABLE_NAME_TO_CLASS_DICT[db_name]
        descriptions = self.session.query(description_class).filter(description_class.id.in_(ids)).all()
        if len(descriptions) == 0:
            warn("No descriptions were found for your id's. Does this %s look like an id from %s" % (list(ids)[0],
                                                                                                     db_name))
        return {i.id: i.description for i in descriptions}

    @staticmethod
    def get_database_names():
        return TABLE_NAME_TO_CLASS_DICT.keys()
