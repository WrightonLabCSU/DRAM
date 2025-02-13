//
// Subworkflow with functionality specific to the WrightonLabCSU/dram pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CALL_GENES                                    } from "${projectDir}/modules/local/call/call_genes_prodigal.nf"
include { QUAST                                         } from "${projectDir}/modules/local/call/quast.nf"
include { QUAST_COLLECT                                 } from "${projectDir}/modules/local/call/quast_collect.nf"
include { GENE_LOCS                                     } from "${projectDir}/modules/local/annotate/gene_locs.nf"

include { GENERIC_HMM_FORMATTER                         } from "${projectDir}/modules/local/annotate/generic_hmm_formatter.nf"  // TODO, This has hard coded paths on the server to the python file. Need to fix this.
include { KEGG_HMM_FORMATTER                            } from "${projectDir}/modules/local/annotate/kegg_hmm_formatter.nf"
include { KOFAM_HMM_FORMATTER                           } from "${projectDir}/modules/local/annotate/kofam_hmm_formatter.nf"
include { DBCAN_HMM_FORMATTER                           } from "${projectDir}/modules/local/annotate/dbcan_hmm_formatter.nf"
include { VOG_HMM_FORMATTER                             } from "${projectDir}/modules/local/annotate/vog_hmm_formatter.nf"
include { CAMPER_HMM_FORMATTER                          } from "${projectDir}/modules/local/annotate/camper_hmm_formatter.nf"
include { CANTHYD_HMM_FORMATTER                         } from "${projectDir}/modules/local/annotate/canthyd_hmm_formatter.nf"
include { SULFUR_HMM_FORMATTER                          } from "${projectDir}/modules/local/annotate/sulfur_hmm_formatter.nf"
include { FEGENIE_HMM_FORMATTER                         } from "${projectDir}/modules/local/annotate/fegenie_hmm_formatter.nf"

include { ADD_TAXA                                      } from "${projectDir}/modules/local/annotate/add_taxa.nf"
include { ADD_BIN_QUALITY                               } from "${projectDir}/modules/local/annotate/add_bin_quality.nf"

include { COMBINE_ANNOTATIONS                           } from "${projectDir}/modules/local/annotate/combine_annotations.nf"
include { COUNT_ANNOTATIONS                             } from "${projectDir}/modules/local/annotate/count_annotations.nf"
include { ADD_ANNOTATIONS                               } from "${projectDir}/modules/local/annotate/add_annotations.nf"
include { MERGE_ANNOTATIONS                             } from "${projectDir}/modules/local/annotate/merge_annotations.nf"
include { GENERATE_GFF_GENBANK                          } from "${projectDir}/modules/local/annotate/generate_gff_genbank.nf"


include { MMSEQS_INDEX                                  } from "${projectDir}/modules/local/annotate/mmseqs_index.nf"

// NextFlow only process with the same name in the same workflow, so either alias it or include it a different workflow
include { MMSEQS_SEARCH as MMSEQS_SEARCH_MEROPS         } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_VIRAL          } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_CAMPER         } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_METHYL         } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_CANTHYD        } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_KEGG           } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_UNIREF         } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_PFAM           } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"

include { ADD_SQL_DESCRIPTIONS as SQL_UNIREF            } from "${projectDir}/modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_VIRAL             } from "${projectDir}/modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_MEROPS            } from "${projectDir}/modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_KEGG              } from "${projectDir}/modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_PFAM              } from "${projectDir}/modules/local/annotate/add_sql_descriptions.nf"

include { HMM_SEARCH as HMM_SEARCH_KOFAM                } from "${projectDir}/modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_DBCAN                } from "${projectDir}/modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_VOG                  } from "${projectDir}/modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_CAMPER               } from "${projectDir}/modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_CANTHYD              } from "${projectDir}/modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_SULFUR               } from "${projectDir}/modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_FEGENIE              } from "${projectDir}/modules/local/annotate/hmmsearch.nf"

include { PARSE_HMM as PARSE_HMM_KOFAM                  } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"
include { PARSE_HMM as PARSE_HMM_DBCAN                  } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"
include { PARSE_HMM as PARSE_HMM_VOG                    } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"
include { PARSE_HMM as PARSE_HMM_CAMPER                 } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"
include { PARSE_HMM as PARSE_HMM_CANTHYD                } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"
include { PARSE_HMM as PARSE_HMM_SULFUR                 } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"
include { PARSE_HMM as PARSE_HMM_FEGENIE                } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO ANNOTATE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ADD_AND_COMBINE {
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

    // Check for additional user-provided annotations
    if( params.add_annotations ){
        ch_add_annots = file(params.add_annotations)
        ADD_ANNOTATIONS( ch_updated_taxa_annots, ch_add_annots )
        ch_final_annots = ADD_ANNOTATIONS.out.combined_annots_out

        // // If the user wants to run trees, do it before we count the annotations
        // if( params.trees ){
        //     TREES( ch_final_annots, params.trees_list, ch_collected_faa, ch_tree_data_files, ch_trees_scripts, ch_add_trees )
        // }

        COUNT_ANNOTATIONS ( ch_final_annots )
        ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts
        ch_annotations_sqlite3 = COUNT_ANNOTATIONS.out.annotations_sqlite3
    }
    else{
        // If the user wants to run trees, do it before we count the annotations
        // if( params.trees ){
        //     TREES( ch_updated_taxa_annots, params.trees_list, ch_collected_faa, ch_tree_data_files, ch_trees_scripts, ch_add_trees )
        // }

        ch_final_annots = ch_updated_taxa_annots
        COUNT_ANNOTATIONS ( ch_final_annots )
        ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts
        ch_annotations_sqlite3 = COUNT_ANNOTATIONS.out.annotations_sqlite3
    }

    if( params.generate_gff || params.generate_gbk ){

        if (!params.call) {
            ch_called_genes = Channel
                .fromPath(file(params.input_genes) / params.genes_fna_fmt, checkIfExists: true)
                .ifEmpty { exit 1, "If you specify --generate_gff or --generate_gbk without --call, you must provide a fasta file of called genes using --input_genes and --genes_fna_fmt,. Cannot find any called gene fasta files matching: ${params.input_genes} and ${params.genes_fna_fmt}\nNB: Path needs to follow pattern: path/to/directory/" }
                .map {
                    sampleName = it.getName().replaceAll(/\.[^.]+$/, '').replaceAll(/\./, '-')
                    tuple(sampleName, it)
                }
            // Collect all individual fasta to pass to quast
            Channel.empty()
                .mix( ch_called_genes  )
                .collect()
                .set { ch_collected_fna }
        }
    
        GENERATE_GFF_GENBANK( ch_collected_fna, params.database_list, ch_final_annots )
    }

}
