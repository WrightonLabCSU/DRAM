//
// Subworkflow with functionality specific to the WrightonLabCSU/dram pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ADD_TAXA                                      } from "${projectDir}/modules/local/add_and_combine/add_taxa.nf"
include { ADD_BIN_QUALITY                               } from "${projectDir}/modules/local/add_and_combine/add_bin_quality.nf"

include { COUNT_ANNOTATIONS                             } from "${projectDir}/modules/local/add_and_combine/count_annotations.nf"
include { ADD_ANNOTATIONS                               } from "${projectDir}/modules/local/add_and_combine/add_annotations.nf"
include { GENERATE_GFF_GENBANK                          } from "${projectDir}/modules/local/add_and_combine/generate_gff_genbank.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO ADD_AND_COUNT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ADD_AND_COUNT {
    take:
    ch_combined_annotations  // channel: [ path(combined_annotations_out) ]

    main:




    // Add Bin Quality to annotations
    if( params.bin_quality ){
        ch_bin_quality = file(params.bin_quality)
        ADD_BIN_QUALITY( ch_combined_annotations, ch_bin_quality )
        ch_updated_annots = ADD_BIN_QUALITY.out.annots_bin_quality_out
    }
    else{
        ch_updated_annots = ch_combined_annotations
    }

    // Add Taxonomy to annotations
    if( params.taxa ){
        ch_taxa = file(params.taxa)
        ADD_TAXA( ch_updated_annots, ch_taxa )
        ch_updated_taxa_annots = ADD_TAXA.out.annots_taxa_out
    }
    else{
        ch_updated_taxa_annots = ch_combined_annotations
    }

    ch_final_annots = ch_updated_taxa_annots
    // Check for additional user-provided annotations
    if( params.add_annotations ){
        ch_add_annots = file(params.add_annotations)
        ADD_ANNOTATIONS( ch_updated_taxa_annots, ch_add_annots )
        ch_final_annots = ADD_ANNOTATIONS.out.combined_annots_out

        // // If the user wants to run trees, do it before we count the annotations
        // if( params.trees ){
        //     TREES( ch_final_annots, params.trees_list, ch_collected_faa, ch_tree_data_files, ch_trees_scripts, ch_add_trees )
        // }


        // Not running count after changing distill to not use sql db built here, remove soon once sure
        // COUNT_ANNOTATIONS ( ch_final_annots )
        // ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts
        // ch_annotations_sqlite3 = COUNT_ANNOTATIONS.out.annotations_sqlite3
    }
    // else{
        // If the user wants to run trees, do it before we count the annotations
        // if( params.trees ){
        //     TREES( ch_updated_taxa_annots, params.trees_list, ch_collected_faa, ch_tree_data_files, ch_trees_scripts, ch_add_trees )
        // }

        // ch_final_annots = ch_updated_taxa_annots
        // COUNT_ANNOTATIONS ( ch_final_annots )
        // ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts
        // ch_annotations_sqlite3 = COUNT_ANNOTATIONS.out.annotations_sqlite3
    // }

    if( params.generate_gff || params.generate_gbk ){

        if (!params.call) {
            ch_called_genes = Channel
                .fromPath(file(params.input_genes) / params.genes_fna_fmt, checkIfExists: true)
                .ifEmpty { exit 1, "If you specify --generate_gff or --generate_gbk without --call, you must provide a fasta file of called genes using --input_genes and --genes_fna_fmt,. Cannot find any called gene fasta files matching: ${params.input_genes} and ${params.genes_fna_fmt}\nNB: Path needs to follow pattern: path/to/directory/" }
                .map {
                    input_fastaName = it.getName().replaceAll(/\.[^.]+$/, '').replaceAll(/\./, '-')
                    tuple(input_fastaName, it)
                }
            // Collect all individual fasta to pass to quast
            Channel.empty()
                .mix( ch_called_genes  )
                .collect()
                .set { ch_collected_fna }
        }
    
        GENERATE_GFF_GENBANK( ch_collected_fna, params.database_list, ch_final_annots )
    }

    emit:
    ch_final_annots

}
