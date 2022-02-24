from os import path, remove
from pkg_resources import resource_filename
import json
import gzip
from shutil import copy2
import warnings

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import pandas as pd

from mag_annotator.database_setup import TABLE_NAME_TO_CLASS_DICT, create_description_db
from mag_annotator.utils import divide_chunks

SEARCH_DATABASES = ('kegg', 'kofam', 'kofam_ko_list', 'uniref', 'pfam', 'dbcan', 'viral', 'peptidase', 'vogdb')
DRAM_SHEETS = ('genome_summary_form', 'module_step_form', 'etc_module_database', 'function_heatmap_form',
               'amg_database')
DATABASE_DESCRIPTIONS = ('pfam_hmm_dat', 'dbcan_fam_activities', 'vog_annotations')

# TODO: store all sequence db locs within database handler class
# TODO: store scoring information here e.g. bitscore_threshold, hmm cutoffs
# TODO: set up custom databases here
# TODO: in advanced config separate search databases, search database description files, description db, DRAM sheets
# TODO: ko_list should be parsed into the DB and stored as a database description file and not a search database


def get_config_loc():
    return path.abspath(resource_filename('mag_annotator', 'CONFIG'))


class DatabaseHandler:
    def __init__(self, config_loc=None):
        if config_loc is None:
            config_loc = get_config_loc()

        # read in configuration # TODO: validate config file after reading it in
        self.config_loc = config_loc
        config = json.loads(open(config_loc).read())
        self.db_locs = {key: value for key, value in config.items() if key in SEARCH_DATABASES}
        self.db_description_locs = {key: value for key, value in config.items() if key in DATABASE_DESCRIPTIONS}
        self.dram_sheet_locs = {key: value for key, value in config.items() if key in DRAM_SHEETS}

        # set up description database connection
        self.description_loc = config.get('description_db')
        if self.description_loc is None:
            self.session = None
            warnings.warn('Database does not exist at path %s' % self.description_loc)
        elif not path.exists(self.description_loc):
            self.session = None
            warnings.warn('Database does not exist at path %s' % self.description_loc)
        else:
            self.start_db_session()

    def start_db_session(self):
        engine = create_engine('sqlite:///%s' % self.description_loc)
        db_session = sessionmaker(bind=engine)
        self.session = db_session()

    # functions for adding descriptions to tables
    def add_descriptions_to_database(self, description_list, db_name, clear_table=True):
        description_class = TABLE_NAME_TO_CLASS_DICT[db_name]
        if clear_table:
            self.session.query(description_class).delete()
        self.session.bulk_save_objects([description_class(**i) for i in description_list])  # TODO: try batching
        self.session.commit()
        self.session.expunge_all()

    # functions for getting descriptions from tables
    def get_description(self, annotation_id, db_name):
        return self.session.query(TABLE_NAME_TO_CLASS_DICT[db_name]).filter_by(id=annotation_id).one().description

    def get_descriptions(self, ids, db_name):
        description_class = TABLE_NAME_TO_CLASS_DICT[db_name]
        descriptions = [
            des 
            for chunk in divide_chunks(list(ids), 499)
            for des in self.session.query(description_class).filter(description_class.id.in_(chunk)).all() 
        ]
        if len(descriptions) == 0:
            warnings.warn("No descriptions were found for your id's. Does this %s look like an id from %s" % (list(ids)[0],
                                                                                                     db_name))
        return {i.id: i.description for i in descriptions}

    @staticmethod
    def get_database_names():
        return TABLE_NAME_TO_CLASS_DICT.keys()

    def set_database_paths(self, kegg_db_loc=None, kofam_hmm_loc=None, kofam_ko_list_loc=None, uniref_db_loc=None,
                           pfam_db_loc=None, pfam_hmm_dat=None, dbcan_db_loc=None, dbcan_fam_activities=None,
                           viral_db_loc=None, peptidase_db_loc=None, vogdb_db_loc=None, vog_annotations=None,
                           description_db_loc=None, genome_summary_form_loc=None, module_step_form_loc=None,
                           etc_module_database_loc=None, function_heatmap_form_loc=None, amg_database_loc=None,
                           write_config=True):
        def check_exists_and_add_to_location_dict(loc, old_value):
            if loc is None:  # if location is none then return the old value
                return old_value
            elif path.isfile(loc):  # if location exists return full path
                return path.realpath(loc)
            else:  # if the location doesn't exist then raise error
                raise ValueError("Database location does not exist: %s" % loc)
        self.db_locs['kegg'] = check_exists_and_add_to_location_dict(kegg_db_loc, self.db_locs.get('kegg'))
        self.db_locs['kofam'] = check_exists_and_add_to_location_dict(kofam_hmm_loc, self.db_locs.get('kofam'))
        self.db_locs['kofam_ko_list'] = check_exists_and_add_to_location_dict(kofam_ko_list_loc,
                                                                              self.db_locs.get('kofam_ko_list'))
        self.db_locs['uniref'] = check_exists_and_add_to_location_dict(uniref_db_loc, self.db_locs.get('uniref'))
        self.db_locs['pfam'] = check_exists_and_add_to_location_dict(pfam_db_loc, self.db_locs.get('pfam'))

        self.db_locs['dbcan'] = check_exists_and_add_to_location_dict(dbcan_db_loc, self.db_locs.get('dbcan'))

        self.db_locs['viral'] = check_exists_and_add_to_location_dict(viral_db_loc, self.db_locs.get('viral'))
        self.db_locs['peptidase'] = check_exists_and_add_to_location_dict(peptidase_db_loc,
                                                                          self.db_locs.get('peptidase'))
        self.db_locs['vogdb'] = check_exists_and_add_to_location_dict(vogdb_db_loc, self.db_locs.get('vogdb'))

        self.db_description_locs['pfam_hmm_dat'] = \
            check_exists_and_add_to_location_dict(pfam_hmm_dat, self.db_description_locs.get('pfam_hmm_dat'))
        self.db_description_locs['dbcan_fam_activities'] = \
            check_exists_and_add_to_location_dict(dbcan_fam_activities,
                                                  self.db_description_locs.get('dbcan_fam_activities'))
        self.db_description_locs['vog_annotations'] = \
            check_exists_and_add_to_location_dict(vog_annotations, self.db_description_locs.get('vog_annotations'))

        self.dram_sheet_locs['genome_summary_form'] = \
            check_exists_and_add_to_location_dict(genome_summary_form_loc,
                                                  self.dram_sheet_locs.get('genome_summary_form'))
        self.dram_sheet_locs['module_step_form'] = \
            check_exists_and_add_to_location_dict(module_step_form_loc, self.dram_sheet_locs.get('module_step_form'))
        self.dram_sheet_locs['etc_module_database'] = \
            check_exists_and_add_to_location_dict(etc_module_database_loc,
                                                  self.dram_sheet_locs.get('etc_module_database'))
        self.dram_sheet_locs['function_heatmap_form'] = \
            check_exists_and_add_to_location_dict(function_heatmap_form_loc,
                                                  self.dram_sheet_locs.get('function_heatmap_form'))
        self.dram_sheet_locs['amg_database'] = \
            check_exists_and_add_to_location_dict(amg_database_loc, self.dram_sheet_locs.get('amg_database'))

        self.description_loc = check_exists_and_add_to_location_dict(description_db_loc, self.description_loc)
        self.start_db_session()

        if write_config:
            self.write_config()

    def write_config(self, config_loc=None):
        if config_loc is None:
            config_loc = self.config_loc
        with open(config_loc, 'w') as f:
            config = {key: value for dict_ in (self.db_locs, self.db_description_locs, self.dram_sheet_locs)
                      for key, value in dict_.items()}
            config['description_db'] = self.description_loc
            f.write(json.dumps(config))

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
            warnings.warn(f'There are {len(mmseqs_ids_not_unique)} non unique headers in {mmseqs_db}! You should definitly investigate this!')
        return mmseqs_headers_split

    @staticmethod
    def process_pfam_descriptions(pfam_hmm_dat):
        if pfam_hmm_dat.endswith('.gz'):
            f = gzip.open(pfam_hmm_dat, 'r').read().decode('utf-8')
        else:
            f = open(pfam_hmm_dat).read()
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
    def process_dbcan_descriptions(dbcan_fam_activities):
        f = open(dbcan_fam_activities)
        description_list = list()
        for line in f.readlines():
            if not line.startswith('#') and len(line.strip()) != 0:
                line = line.strip().split()
                if len(line) == 1:
                    description = line[0]
                elif line[0] == line[1]:
                    description = ' '.join(line[1:])
                else:
                    description = ' '.join(line)
                description_list.append({'id': line[0], 'description': description.replace('\n', ' ')})
        return description_list

    @staticmethod
    def process_vogdb_descriptions(vog_annotations):
        annotations_table = pd.read_csv(vog_annotations, sep='\t', index_col=0)
        annotations_list = [{'id': vog, 'description': '%s; %s' % (row['ConsensusFunctionalDescription'],
                                                                   row['FunctionalCategory'])}
                            for vog, row in annotations_table.iterrows()]
        return annotations_list

    # TODO: Make option to build on description database that already exists?
    def populate_description_db(self, output_loc=None, update_config=True):
        if self.description_loc is None and output_loc is None:  # description db location must be set somewhere
            raise ValueError('Must provide output location if description db location is not set in configuration')
        if output_loc is not None:  # if new description db location is set then save it there
            self.description_loc = output_loc
            self.start_db_session()
        if path.exists(self.description_loc):
            remove(self.description_loc)
        create_description_db(self.description_loc)

        # fill database
        if self.db_locs.get('kegg') is not None:
            self.add_descriptions_to_database(self.make_header_dict_from_mmseqs_db(self.db_locs['kegg']), 'kegg_description',
                                              clear_table=True)
        if self.db_locs.get('uniref') is not None:
            self.add_descriptions_to_database(self.make_header_dict_from_mmseqs_db(self.db_locs['uniref']) ,
                                              'uniref_description', clear_table=True)
        if self.db_description_locs.get('pfam_hmm_dat') is not None:
            self.add_descriptions_to_database(self.process_pfam_descriptions(self.db_description_locs['pfam_hmm_dat']),
                                              'pfam_description', clear_table=True)
        if self.db_description_locs.get('dbcan_fam_activities') is not None:
            self.add_descriptions_to_database(self.process_dbcan_descriptions(
                self.db_description_locs['dbcan_fam_activities']), 'dbcan_description', clear_table=True)
        if self.db_locs.get('viral') is not None:
            self.add_descriptions_to_database(self.make_header_dict_from_mmseqs_db(self.db_locs['viral']),
                                              'viral_description', clear_table=True)
        if self.db_locs.get('peptidase') is not None:
            self.add_descriptions_to_database(self.make_header_dict_from_mmseqs_db(self.db_locs['peptidase']),
                                              'peptidase_description', clear_table=True)
        if self.db_description_locs.get('vog_annotations') is not None:
            self.add_descriptions_to_database(
                self.process_vogdb_descriptions(self.db_description_locs['vog_annotations']), 'vogdb_description',
                clear_table=True)

        if update_config:  # if new description db is set then save it
            self.write_config()

    def print_database_locations(self):
        # search databases
        print('Processed search databases')
        print('KEGG db: %s' % self.db_locs.get('kegg'))
        print('KOfam db: %s' % self.db_locs.get('kofam'))
        print('KOfam KO list: %s' % self.db_locs.get('kofam_ko_list'))
        print('UniRef db: %s' % self.db_locs.get('uniref'))
        print('Pfam db: %s' % self.db_locs.get('pfam'))
        print('dbCAN db: %s' % self.db_locs.get('dbcan'))
        print('RefSeq Viral db: %s' % self.db_locs.get('viral'))
        print('MEROPS peptidase db: %s' % self.db_locs.get('peptidase'))
        print('VOGDB db: %s' % self.db_locs.get('vogdb'))
        print()
        # database descriptions used during description db population
        print('Descriptions of search database entries')
        print('Pfam hmm dat: %s' % self.db_description_locs.get('pfam_hmm_dat'))
        print('dbCAN family activities: %s' % self.db_description_locs.get('dbcan_fam_activities'))
        print('VOG annotations: %s' % self.db_description_locs.get('vog_annotations'))
        print()
        # description database
        print('Description db: %s' % self.description_loc)
        print()
        # DRAM sheets
        print('DRAM distillation sheets')
        print('Genome summary form: %s' % self.dram_sheet_locs.get('genome_summary_form'))
        print('Module step form: %s' % self.dram_sheet_locs.get('module_step_form'))
        print('ETC module database: %s' % self.dram_sheet_locs.get('etc_module_database'))
        print('Function heatmap form: %s' % self.dram_sheet_locs.get('function_heatmap_form'))
        print('AMG database: %s' % self.dram_sheet_locs.get('amg_database'))

    def filter_db_locs(self, low_mem_mode=False, use_uniref=True, use_vogdb=True, master_list=None):
        if master_list is None:
            dbs_to_use = self.db_locs.keys()
        else:
            dbs_to_use = master_list
        # filter out dbs for low mem mode
        if low_mem_mode:
            if ('kofam' not in self.db_locs) or ('kofam_ko_list' not in self.db_locs):
                raise ValueError('To run in low memory mode KOfam must be configured for use in DRAM')
            dbs_to_use = [i for i in dbs_to_use if i not in ('uniref', 'kegg', 'vogdb')]
        # check on uniref status
        if use_uniref:
            if 'uniref' not in self.db_locs:
                warnings.warn('Sequences will not be annoated against uniref as it is not configured for use in DRAM')
        else:
            dbs_to_use = [i for i in dbs_to_use if i != 'uniref']
        # check on vogdb status
        if use_vogdb:
            if 'vogdb' not in self.db_locs:
                warnings.warn('Sequences will not be annoated against VOGDB as it is not configured for use in DRAM')
        else:
            dbs_to_use = [i for i in dbs_to_use if i != 'vogdb']
        self.db_locs = {key: value for key, value in self.db_locs.items() if key in dbs_to_use}

    def clear_config(self):
        self.db_locs = {}
        self.db_description_locs = {}
        self.dram_sheet_locs = {}
        self.description_loc = None


def set_database_paths(kegg_db_loc=None, kofam_hmm_loc=None, kofam_ko_list_loc=None, uniref_db_loc=None,
                       pfam_db_loc=None, pfam_hmm_dat=None, dbcan_db_loc=None, dbcan_fam_activities=None,
                       viral_db_loc=None, peptidase_db_loc=None, vogdb_db_loc=None, vog_annotations=None,
                       description_db_loc=None, genome_summary_form_loc=None, module_step_form_loc=None,
                       etc_module_database_loc=None, function_heatmap_form_loc=None, amg_database_loc=None,
                       clear_config=False, update_description_db=False):
    db_handler = DatabaseHandler()
    if clear_config:
        db_handler.clear_config()
    db_handler.set_database_paths(kegg_db_loc, kofam_hmm_loc, kofam_ko_list_loc, uniref_db_loc,
                                  pfam_db_loc, pfam_hmm_dat, dbcan_db_loc, dbcan_fam_activities,
                                  viral_db_loc, peptidase_db_loc, vogdb_db_loc, vog_annotations,
                                  description_db_loc, genome_summary_form_loc, module_step_form_loc,
                                  etc_module_database_loc, function_heatmap_form_loc, amg_database_loc,
                                  write_config=True)
    if update_description_db:
        db_handler.populate_description_db()


def print_database_locations(config_loc=None):
    db_handler = DatabaseHandler(config_loc)
    db_handler.print_database_locations()


def populate_description_db(output_loc=None, config_loc=None):
    db_handler = DatabaseHandler(config_loc)
    db_handler.populate_description_db(output_loc)


def export_config(output_file=None):
    config_loc = get_config_loc()
    if output_file is None:
        print(open(config_loc).read())
    else:
        copy2(config_loc, output_file)


def import_config(config_loc):
    system_config = get_config_loc()
    copy2(config_loc, system_config)
