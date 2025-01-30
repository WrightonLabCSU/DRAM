//
// Subworkflow with functionality specific to the WrightonLabCSU/dram pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MERGE_ANNOTATIONS                                    } from "${projectDir}/modules/local/annotate/merge_annotations.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO MERGE ANNOTATIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MERGE {

    // Verify the directory exists
    def annotations_dir = file(params.merge_annotations)
    if (!annotations_dir.exists()) {
        error "Error: The specified directory for merging annotations (--merge_annotations) does not exist: ${params.merge_annotations}"
    }

    // Verify the directory contains .tsv files
    def tsv_files = annotations_dir.list().findAll { it.endsWith('.tsv') }
    if (tsv_files.isEmpty()) {
        error "Error: The specified directory for merging annotations (--merge_annotations) does not contain any .tsv files: ${params.merge_annotations}"
    }

    // Create a channel with the paths to the .tsv files
    Channel
        .from(tsv_files.collect { annotations_dir.toString() + '/' + it })
        .set { ch_merge_annotations }
    Channel.empty()
        .mix( ch_merge_annotations )
        .collect()
        .set { ch_merge_annotations_collected }

    MERGE_ANNOTATIONS( ch_merge_annotations_collected )
}