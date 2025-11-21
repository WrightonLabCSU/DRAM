include { logColours           } from '../nf-core/utils_nfcore_pipeline'
include { getWorkflowVersion           } from '../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


def getDBFlag(db_list, db_name, value_for_all) {
    if (db_list.contains(value_for_all)) {
        return true
    } else if (db_list.contains(db_name)) {
        return true
    } else {
        return false
    }
}