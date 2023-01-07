"""
One test to stop so much pain, did you accidentally commit a test config?
This tests that the config is clean. It will fail when you are testing a working copy but not when you push to github.

"""

import pytest
import logging
from os import path
from mag_annotator.database_handler import DatabaseHandler
from mag_annotator.utils import setup_logger



def test_clean_config():
    empty_config = {
        "search_databases": {
        "kegg": None,
        "kofam_hmm": None,
        "kofam_ko_list": None,
        "uniref": None,
        "pfam": None,
        "dbcan": None,
        "viral": None,
        "peptidase": None,
        "vogdb": None
      },
      "custom_dbs": None,
      "database_descriptions": {
        "pfam_hmm_dat": None,
        "dbcan_fam_activities": None,
        "vog_annotations": None
      },
      "dram_sheets": {
        "genome_summary_form": None,
        "module_step_form": None,
        "etc_module_database": None,
        "function_heatmap_form": None,
        "amg_database": None
      },
      "description_db": None,
      "dram_version": None
    }
    logger = logging.getLogger('test_log')
    db_handler = DatabaseHandler(logger)
    assert db_handler.config == empty_config
