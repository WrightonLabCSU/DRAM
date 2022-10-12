from os import path, remove, getenv
from pkg_resources import resource_filename
import json
import gzip
import logging
from shutil import copy2
import warnings
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from datetime import datetime
from functools import partial

import pandas as pd

from mag_annotator import __version__ as current_dram_version
from mag_annotator.database_setup import TABLE_NAME_TO_CLASS_DICT, create_description_db
from mag_annotator.utils import divide_chunks, setup_logger

SEARCH_DATABASES = {'kegg', 'kofam_hmm', 'kofam_ko_list', 'uniref', 'pfam', 'dbcan',
                    'viral', 'peptidase', 'vogdb' }
DRAM_SHEETS = ('genome_summary_form', 'module_step_form', 'etc_module_database', 'function_heatmap_form',
                'amg_database')
DATABASE_DESCRIPTIONS = ('pfam_hmm', 'dbcan_fam_activities', 'vog_annotations')

# TODO: store all sequence db locs within database handler class
# TODO: store scoring information here e.g. bitscore_threshold, hmm cutoffs
# TODO: set up custom databases here
# TODO: in advanced config separate search databases, search database description files, description db, DRAM sheets
# TODO: ko_list should be parsed into the DB and stored as a database description file and not a search database


def get_config_loc():
    loc = getenv('DRAM_CONFIG_LOCATION')
    if loc:
        return loc
    else:
        return path.abspath(resource_filename('mag_annotator', 'CONFIG'))


def clear_dict(val):
    if isinstance(val, dict):
        return {k: clear_dict(v) for k, v in val.items()}
    else:
        return None


class DatabaseHandler:

    def __init__(self, logger, config_loc=None):
        # read in new configuration
        # TODO: validate config file after reading it in
        if logger is None:
            logger = logging.getLogger("database_handler.log")
            # log_path = self.get_log_path()
            # setup_logger(logger, log_path)
            setup_logger(logger)
            logger.info(f"Logging to console")
        self.logger = logger
        if config_loc is None:
            config_loc = get_config_loc()
        self.load_config(config_loc)

        self.config_loc = config_loc



    def load_config(self, config_file):
        conf = json.loads(open(config_file).read())
        if len(conf) == 0:
            self.logger.warn('There is no config information in the provided file')
            self.clear_config(write_config=False)
        if 'dram_version' not in conf:
            warnings.warn("The DRAM version in your config is empty."
                          " This may not be a problem, but if this"
                          " import fails then you should check that"
                          " the origin of the file is valid.")
            self.__construct_from_dram_pre_1_4_0(conf)
        else:
            conf_version = conf.get('dram_version')
            if conf_version is None:
                self.__construct_from_dram_pre_1_4_0(conf)
            elif conf_version not in {current_dram_version, "1.4.0",  "1.4.0rc1", "1.4.0rc2", "1.4.0rc3", "1.4.0rc4"}: # Known suported versions
                warnings.warn("The DRAM version in your config is not listed in the versions "
                              "that are known to work. This may not be a problem, but if this "
                              "import fails then you should contact suport.")
            self.__construct_default(conf)

    def get_log_path(self):
        path = self.config.get('log_path')
        if path is None:
            path = path.join(self.config_loc, 'database_processing.log')
        return path


    def __construct_default(self, conf:dict):
        self.config = conf

        # set up description database connection
        description_loc = self.config.get('description_db')
        if description_loc is None:
            self.session = None
            warnings.warn('Database does not exist at path %s' % description_loc)
        elif not path.exists(description_loc):
            self.session = None
            warnings.warn('Database does not exist at path %s' % description_loc)
        else:
            self.start_db_session()

    def __construct_from_dram_pre_1_4_0(self, config_old):
        """
        Import older dram configs that predate 1.3


        :param config_old: A config with no dram version so older than 1.3
        """

        system_config_loc = get_config_loc()
        self.config = clear_dict(json.loads(open(system_config_loc).read()))
        self.config_loc = system_config_loc


        # read in configuration # TODO: validate config file after reading it in
        self.config['search_databases'] = {
            key: value for key, value in config_old.items() if key in SEARCH_DATABASES}
        self.config['database_descriptions']  = {
            key: value for key, value in config_old.items() if key in DATABASE_DESCRIPTIONS}
        self.config['dram_sheets']  = {
            key: value for key, value in config_old.items() if key in DRAM_SHEETS}
        self.config["dram_version"] = current_dram_version

        # set up description database connection
        self.config['description_db'] = config_old.get('description_db')
        if self.config.get('description_db') is None:
            self.session = None
            warnings.warn('Database does not exist at path %s' % self.config.get('description_db'))
        elif not path.exists(self.config.get('description_db')):
            self.session = None
            warnings.warn('Database does not exist at path %s' % self.config.get('description_db'))
        else:
            self.start_db_session()

    def start_db_session(self):
        engine = create_engine('sqlite:///%s' % self.config.get('description_db'))
        db_session = sessionmaker(bind=engine)
        self.session = db_session()


    # functions for adding descriptions to tables
    def add_descriptions_to_database(self, description_list, db_name, clear_table=True):
        description_class = TABLE_NAME_TO_CLASS_DICT[db_name]
        if clear_table:
            self.session.query(description_class).delete()
        # TODO: try batching
        self.session.bulk_save_objects([description_class(**i) for i in description_list])
        self.session.commit()
        self.session.expunge_all()

    # functions for getting descriptions from tables
    def get_description(self, annotation_id, db_name, return_ob=False):
        return self.session.query(TABLE_NAME_TO_CLASS_DICT[db_name]).filter_by(id=annotation_id).one().description

    def get_descriptions(self, ids, db_name, description_name='description'):
        description_class = TABLE_NAME_TO_CLASS_DICT[db_name]
        descriptions = [
            des
            for chunk in divide_chunks(list(ids), 499)
            for des in self.session.query(description_class).filter(description_class.id.in_(chunk)).all()
        ]
        # [des for des in self.session.query(description_class).filter(description_class.id.in_(list(ids))).all() ]
        # [i.id for i in self.session.query(TABLE_NAME_TO_CLASS_DICT['dbcan_description']).all()]
        if len(descriptions) == 0:
            warnings.warn("No descriptions were found for your id's. Does this %s look like an id from %s" % (list(ids)[0],
                                                                                                     db_name))
        return {i.id: i.__dict__[description_name] for i in descriptions}

    @staticmethod
    def get_database_names():
        return TABLE_NAME_TO_CLASS_DICT.keys()

    def get_settings_str(self):
        out_str = ""
        settings = self.config.get('setup_info')
        if settings is None:
            warnings.warn('there are no settings, the config is corrupted or too old.', DeprecationWarning)
            return 'there are no settings, the config is corrupted or too old.'
        for i in ["search_databases", "database_descriptions", "dram_sheets"]:
            out_str += "\n"
            for k in self.config.get(i):
                if settings.get(k) is not None:
                    out_str += f"\n{settings[k]['name']}:"
                    for l, w in settings[k].items():
                        if l =='name':
                            continue
                        out_str += f"\n    {l.title()}: {w}"
        return out_str


    def set_database_paths(self, kegg_loc=None, kofam_hmm_loc=None, kofam_ko_list_loc=None, uniref_loc=None,
                           pfam_loc=None, pfam_hmm_loc=None, dbcan_loc=None, dbcan_fam_activities_loc=None,
                           dbcan_subfam_ec_loc=None, viral_loc=None, peptidase_loc=None, vogdb_loc=None,
                           vog_annotations_loc=None, description_db_loc=None, log_path_loc=None,
                           genome_summary_form_loc=None,
                           module_step_form_loc=None, etc_module_database_loc=None,
                           function_heatmap_form_loc=None, amg_database_loc=None, write_config=True):
        def check_exists_and_add_to_location_dict(loc, old_value):
            if loc is None:  # if location is none then return the old value
                return old_value
            if path.isfile(loc):  # if location exists return full path
                return path.realpath(loc)
            else:  # if the location doesn't exist then raise error
                raise ValueError("Database location does not exist: %s" % loc)
        locs = {
                "search_databases": {
                  'kegg': kegg_loc,
                  'kofam_hmm': kofam_hmm_loc,
                  'kofam_ko_list': kofam_ko_list_loc,
                  'uniref': uniref_loc,
                  'pfam': pfam_loc,
                  'dbcan': dbcan_loc,
                  'viral': viral_loc,
                  'peptidase': peptidase_loc,
                  'vogdb': vogdb_loc,
                },
                "database_descriptions": {
                  'pfam_hmm': pfam_hmm_loc,
                  'dbcan_fam_activities': dbcan_fam_activities_loc,
                  'dbcan_subfam_ec': dbcan_subfam_ec_loc,
                  'vog_annotations': vog_annotations_loc,
                },
                "dram_sheets": {
                  'genome_summary_form': genome_summary_form_loc,
                  'module_step_form': module_step_form_loc,
                  'etc_module_database': etc_module_database_loc,
                  'function_heatmap_form': function_heatmap_form_loc,
                  'amg_database': amg_database_loc,
                },
        }

        self.config.update({i:{
            k:check_exists_and_add_to_location_dict(locs[i][k], self.config.get(i).get(k))
            for k in locs[i]} for i in locs})

        self.config['description_db'] = check_exists_and_add_to_location_dict(description_db_loc, self.config.get('description_db'))
        self.config['log_path'] = check_exists_and_add_to_location_dict(log_path_loc, self.config.get('log_path_db'))
        self.start_db_session()

        if write_config:
            self.write_config()

    def write_config(self, config_loc=None):
        if config_loc is None:
            config_loc = self.config_loc
        with open(config_loc, 'w') as f:
            f.write(json.dumps(self.config, indent=2))

    @staticmethod
    def make_header_dict_from_mmseqs_db(mmseqs_db):
        mmseqs_headers_handle = open('%s_h' % mmseqs_db, 'rb')
        mmseqs_headers = mmseqs_headers_handle.read().decode(errors='ignore')
        mmseqs_headers = [i.strip() for i in mmseqs_headers.strip().split('\n\x00') if len(i) > 0]
        mmseqs_headers_split = []
        mmseqs_ids_unique = set()
        mmseqs_ids_not_unique = set()
        # TODO this could be faster with numpy
        for i in mmseqs_headers:
            header = {'id': i.split(' ')[0], 'description': i}
            if header['id'] not in mmseqs_ids_unique:
                mmseqs_headers_split += [header]
                mmseqs_ids_unique.add(header['id'])
            else:
                mmseqs_ids_not_unique.add(header['id'])
        if len(mmseqs_ids_not_unique) > 0:
            warnings.warn(f'There are {len(mmseqs_ids_not_unique)} non unique headers '
                          f'in {mmseqs_db}! You should definitly investigate this!')
        return mmseqs_headers_split

    @staticmethod
    def process_pfam_descriptions(pfam_hmm):
        if pfam_hmm.endswith('.gz'):
            f = gzip.open(pfam_hmm, 'r').read().decode('utf-8')
        else:
            f = open(pfam_hmm).read()
        entries = f.strip().split('//')
        description_list = list()
        for i, entry in enumerate(entries):
            if len(entry) > 0:
                entry = entry.split('\n')
                ascession = None
                description = None
                for line in entry:
                    line = line.strip()
                    if line.startswith('#=GF AC'):
                        ascession = line.split('   ')[-1]
                    if line.startswith('#=GF DE'):
                        description = line.split('   ')[-1]
                description_list.append({'id': ascession, 'description': description})
        return description_list

    @staticmethod
    def process_dbcan_descriptions(dbcan_fam_activities, dbcan_subfam_ec):
        def line_reader(line):
            if not line.startswith('#') and len(line.strip()) != 0:
                line = line.strip().split()
                if len(line) == 1:
                    description = line[0]
                elif line[0] == line[1]:
                    description = ' '.join(line[1:])
                else:
                    description = ' '.join(line)
                return pd.DataFrame({'id': line[0], 'description': description.replace('\n', ' ')}, index=[0])
        with open(dbcan_fam_activities) as f:
            description_data = pd.concat([line_reader(line) for line in f.readlines()])

        ec_data = (pd.read_csv(dbcan_subfam_ec, sep='\t',names=['id', 'id2','ec'], comment='#')[['id', 'ec']]
                           .drop_duplicates())
        ec_data = (pd.concat([ec_data['id'],
                                     ec_data['ec'].str.split('|', expand=True)]
                                    ,axis=1)
                           .melt(id_vars='id',value_name='ec')
                           .dropna(subset=['ec'])[['id', 'ec']]
                           .groupby('id')
                           .apply(lambda x: ','.join(x['ec'].unique()))
                           )
        ec_data = pd.DataFrame(ec_data, columns=['ec']).reset_index()
        data = pd.merge(description_data, ec_data, how='outer', on='id').fillna('')
        return [i.to_dict() for _, i in data.iterrows()]


    @staticmethod
    def process_vogdb_descriptions(vog_annotations):
        annotations_table = pd.read_csv(vog_annotations, sep='\t', index_col=0)
        annotations_list = [{'id': vog, 'description': '%s; %s' % (row['ConsensusFunctionalDescription'],
                                                                   row['FunctionalCategory'])}
                            for vog, row in annotations_table.iterrows()]
        return annotations_list

    # TODO: Make option to build on description database that already exists?
    def populate_description_db(self, output_loc=None, select_db=None, update_config=True,  erase_old_db=False):
        if self.config.get('description_db') is None and output_loc is None:  # description db location must be set somewhere
            self.logger.critical('Must provide output location if description db location is not set in configuration')
            raise ValueError('Must provide output location if description db location is not set in configuration')
        if output_loc is not None:  # if new description db location is set then save it there
            self.config['description_db']= output_loc
            self.start_db_session()
        # I don't think this is needed
        if path.exists(self.config.get('description_db')) and erase_old_db:
            remove(self.config.get('description_db'))
        create_description_db(self.config.get('description_db'))
        def check_db(db_name, db_function):
            # TODO add these sorts of checks to a separate function
            # if self.config.get('search_databases').get(db_name) is None:
            #     return
            # if not path.exists(self.config['search_databases'][db_name]):
            #     logger.warn(f"There is a path for the {db_name} db in the config, but there"
            #                 " is no file at that path. The path is:"
            #                 f"{self.config['search_databases'][db_name]}")
            #     return
            self.add_descriptions_to_database(
                db_function(),
                f'{db_name}_description',
                clear_table=True)
            self.config['setup_info'][db_name]['description_db_updated'] = \
                datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
            self.logger.info(f'Description updated for the {db_name} database')
        # fill database
        mmseqs_database = ['kegg', 'uniref',  'viral', 'peptidase']
        process_functions = {i:partial(self.make_header_dict_from_mmseqs_db,
                                       self.config['search_databases'][i])
                             for i in mmseqs_database
                             if self.config['search_databases'][i] is not None}
        # Use table names
        process_functions.update({
            'pfam': partial(self.process_pfam_descriptions,
                                self.config.get('database_descriptions')['pfam_hmm']),
            'dbcan': partial(self.process_dbcan_descriptions,
                                            self.config.get('database_descriptions')['dbcan_fam_activities'],
                                            self.config.get('database_descriptions')['dbcan_subfam_ec']),
            'vogdb': partial(self.process_vogdb_descriptions,
                                         self.config.get('database_descriptions')['vog_annotations'])
        })
        if select_db is not None:
            process_functions = {i:k for i, k in process_functions.items() if i in select_db}

        for i, k in process_functions.items():
            check_db(i, k)

        if update_config:  # if new description db is set then save it
            self.write_config()

    def filter_db_locs(self, low_mem_mode=False, use_uniref=True, use_vogdb=True, master_list=None):
        if master_list is None:
            dbs_to_use = self.config['search_databases'].keys()
        else:
            dbs_to_use = master_list
        # filter out dbs for low mem mode
        if low_mem_mode:
            if ('kofam_hmm' not in self.config.get('search_databases')) or ('kofam_ko_list' not in self.config.get('search_databases')):
                raise ValueError('To run in low memory mode KOfam must be configured for use in DRAM')
            dbs_to_use = [i for i in dbs_to_use if i not in ('uniref', 'kegg', 'vogdb')]
        # check on uniref status
        if use_uniref:
            if 'uniref' not in self.config.get('search_databases'):
                warnings.warn('Sequences will not be annoated against uniref as it is not configured for use in DRAM')
        else:
            dbs_to_use = [i for i in dbs_to_use if i != 'uniref']
        # check on vogdb status
        if use_vogdb:
            if 'vogdb' not in self.config.get('search_databases'):
                warnings.warn('Sequences will not be annoated against VOGDB as it is not configured for use in DRAM')
        else:
            dbs_to_use = [i for i in dbs_to_use if i != 'vogdb']
        self.config['search_databases'] = {key: value for key, value in self.config.get('search_databases').items() if key in dbs_to_use}

    def clear_config(self, write_config=False):
        self.config = {
                        "search_databases": {},
                        "database_descriptions": {},
                        "dram_sheets": {},
                        "dram_version": current_dram_version,
                        "description_db": None,
                        "setup_info": {},
                        "log_path": None
                      }
        if write_config:
            self.write_config()


def set_database_paths(clear_config=False, update_description_db=False, **kargs):
    #TODO Add tests
    db_handler = DatabaseHandler(None)
    if clear_config:
        db_handler.clear_config(write_config=True)
    db_handler.set_database_paths(**kargs, write_config=True)
    if update_description_db:
        db_handler.populate_description_db()


def print_database_locations(config_loc=None):
    conf = DatabaseHandler(None, config_loc)
    # search databases
    print('Processed search databases')
    print('KEGG db: %s' % conf.config.get('search_databases').get('kegg'))
    print('KOfam db: %s' % conf.config.get('search_databases').get('kofam_hmm'))
    print('KOfam KO list: %s' % conf.config.get('search_databases').get('kofam_ko_list'))
    print('UniRef db: %s' % conf.config.get('search_databases').get('uniref'))
    print('Pfam db: %s' % conf.config.get('search_databases').get('pfam'))
    print('dbCAN db: %s' % conf.config.get('search_databases').get('dbcan'))
    print('RefSeq Viral db: %s' % conf.config.get('search_databases').get('viral'))
    print('MEROPS peptidase db: %s' % conf.config.get('search_databases').get('peptidase'))
    print('VOGDB db: %s' % conf.config.get('search_databases').get('vogdb'))
    # database descriptions used during description db population
    print('Descriptions of search database entries')
    print('Pfam hmm dat: %s' % conf.config.get('database_descriptions').get('pfam_hmm'))
    print('dbCAN family activities: %s' % conf.config.get('database_descriptions').get('dbcan_fam_activities'))
    print('VOG annotations: %s' % conf.config.get('database_descriptions').get('vog_annotations'))
    print()
    # description database
    print('Description db: %s' % conf.config.get('description_db'))
    print()
    # DRAM sheets
    print('DRAM distillation sheets')
    print('Genome summary form: %s' % conf.config.get('dram_sheets').get('genome_summary_form'))
    print('Module step form: %s' % conf.config.get('dram_sheets').get('module_step_form'))
    print('ETC module database: %s' % conf.config.get('dram_sheets').get('etc_module_database'))
    print('Function heatmap form: %s' % conf.config.get('dram_sheets').get('function_heatmap_form'))
    print('AMG database: %s' % conf.config.get('dram_sheets').get('amg_database'))


def print_database_settings(config_loc=None):
    conf = DatabaseHandler(None, config_loc)
    print(conf.get_settings_str())



def populate_description_db(output_loc=None, select_db=None,  config_loc=None):
    db_handler = DatabaseHandler(None, config_loc)
    db_handler.populate_description_db(output_loc, select_db)


def export_config(output_file=None):
    config_loc = get_config_loc()
    if output_file is None:
        print(open(config_loc).read())
    else:
        copy2(config_loc, output_file)


def import_config(config_loc):
    system_config = get_config_loc()
    db_handler = DatabaseHandler(None, config_loc)
    with open(system_config, "w") as outfile:
        json.dump(db_handler.config, outfile, indent=2)
    print('Import, appears to be successfull.')


def mv_db_folder(new_location:str='./', old_config_file:str=None):
    new_location = path.abspath(new_location)
    old_config_file = path.abspath(old_config_file)
    db_handler = DatabaseHandler(None)
    if old_config_file is not None:
        db_handler.load_config(old_config_file)
    paths = ["search_databases", "dram_sheets", "database_descriptions"]
    def auto_move_path(k:str, v:str):
        if v is None:
            db_handler.logger.warn(f"The path for {k} was not set, so can't update.")
            return
        new_path = path.join(new_location, path.basename(v))
        if not path.exists(new_path):
            db_handler.logger.warn(f"There is no file at path {new_path},"
                                   f" so no new location will be set for {k}.")
            return
        db_handler.logger.info(f"Moving {k} to {new_path}")
        db_handler.set_database_paths(**{f"{k}_loc": new_path}, write_config=True)
    auto_move_path('log_path', db_handler.config.get('log_path'))
    auto_move_path('description_db', db_handler.config.get('description_db'))
    for i in paths:
        for k, v in db_handler.config.get(i).items():
            auto_move_path(k, v)


