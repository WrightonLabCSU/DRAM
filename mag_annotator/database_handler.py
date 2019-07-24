from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from mag_annotator.database_setup import TABLE_NAME_TO_CLASS_DICT
from mag_annotator.utils import get_database_locs


class DatabaseHandler:
    def __init__(self):
        # get configuration
        config = get_database_locs()
        # check if db exists and not then create
        engine = create_engine('sqlite:///%s' % config['description_db_loc'])
        DBSession = sessionmaker(bind=engine)
        self.session = DBSession()

    def get_description(self, id, db_name):
        return self.session.query(TABLE_NAME_TO_CLASS_DICT[db_name]).filter_by(id=id).one().description

    def get_descriptions(self, ids, db_name):
        description_class = TABLE_NAME_TO_CLASS_DICT[db_name]
        descriptions = self.session.query(description_class).filter(description_class.id.in_(ids)).all()
        return {i.id: i.description for i in descriptions}
