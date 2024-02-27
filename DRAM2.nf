/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DRAM2: <Description>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Homepage of project 
    homePage = <GitHub homepage>
    
    Author of DRAM2 Nextflow pipeline
    author = Reed Woyda, Rory Flynn
    institution = Colorado State University - Wrighton Lab

    Description of project
    description = <Description>

    Main pipeline script
    mainScript = DRAM2.nf
    
    version
    v2.0.1
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Enable Nextflow DSL2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Load Modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { RENAME_FASTA                                  } from './modules/call/rename_fasta.nf'
include { CALL_GENES                                    } from './modules/call/call_genes_prodigal.nf'

include { TRNA_SCAN                                     } from './modules/annotate/trna_scan.nf'
include { RRNA_SCAN                                     } from './modules/annotate/rrna_scan.nf'
include { TRNA_COLLECT                                  } from './modules/annotate/trna_collect.nf'
include { RRNA_COLLECT                                  } from './modules/annotate/rrna_collect.nf'
include { ADD_TAXA                                      } from './modules/annotate/add_taxa.nf'
include { ADD_BIN_QUALITY                               } from './modules/annotate/add_bin_quality.nf'

include { MMSEQS_INDEX                                  } from './modules/annotate/mmseqs_index.nf'
include { MMSEQS_SEARCH as MMSEQS_SEARCH_MEROPS         } from './modules/annotate/mmseqs_search.nf'
include { MMSEQS_SEARCH as MMSEQS_SEARCH_VIRAL          } from './modules/annotate/mmseqs_search.nf'
include { MMSEQS_SEARCH as MMSEQS_SEARCH_CAMPER         } from './modules/annotate/mmseqs_search.nf'
include { MMSEQS_SEARCH as MMSEQS_SEARCH_METHYL         } from './modules/annotate/mmseqs_search.nf'
include { MMSEQS_SEARCH as MMSEQS_SEARCH_CANTHYD        } from './modules/annotate/mmseqs_search.nf'
include { MMSEQS_SEARCH as MMSEQS_SEARCH_KEGG           } from './modules/annotate/mmseqs_search.nf'
include { MMSEQS_SEARCH as MMSEQS_SEARCH_UNIREF         } from './modules/annotate/mmseqs_search.nf'
include { MMSEQS_SEARCH as MMSEQS_SEARCH_PFAM           } from './modules/annotate/mmseqs_search.nf'

include { ADD_SQL_DESCRIPTIONS as SQL_UNIREF            } from './modules/annotate/add_sql_descriptions.nf'
include { ADD_SQL_DESCRIPTIONS as SQL_VIRAL             } from './modules/annotate/add_sql_descriptions.nf'
include { ADD_SQL_DESCRIPTIONS as SQL_MEROPS            } from './modules/annotate/add_sql_descriptions.nf'
include { ADD_SQL_DESCRIPTIONS as SQL_KEGG              } from './modules/annotate/add_sql_descriptions.nf'

include { HMM_SEARCH as HMM_SEARCH_KOFAM                } from './modules/annotate/hmmsearch.nf'
include { PARSE_HMM as PARSE_HMM_KOFAM                  } from './modules/annotate/parse_hmmsearch.nf'

include { HMM_SEARCH as HMM_SEARCH_DBCAN                } from './modules/annotate/hmmsearch.nf'
include { PARSE_HMM as PARSE_HMM_DBCAN                  } from './modules/annotate/parse_hmmsearch.nf'

include { HMM_SEARCH as HMM_SEARCH_VOG                  } from './modules/annotate/hmmsearch.nf'
include { PARSE_HMM as PARSE_HMM_VOG                    } from './modules/annotate/parse_hmmsearch.nf'

include { HMM_SEARCH as HMM_SEARCH_CAMPER               } from './modules/annotate/hmmsearch.nf'
include { PARSE_HMM as PARSE_HMM_CAMPER                 } from './modules/annotate/parse_hmmsearch.nf'

include { HMM_SEARCH as HMM_SEARCH_CANTHYD              } from './modules/annotate/hmmsearch.nf'
include { PARSE_HMM as PARSE_HMM_CANTHYD                } from './modules/annotate/parse_hmmsearch.nf'

include { HMM_SEARCH as HMM_SEARCH_SULFUR               } from './modules/annotate/hmmsearch.nf'
include { PARSE_HMM as PARSE_HMM_SULFUR                 } from './modules/annotate/parse_hmmsearch.nf'

include { HMM_SEARCH as HMM_SEARCH_FEGENIE              } from './modules/annotate/hmmsearch.nf'
include { PARSE_HMM as PARSE_HMM_FEGENIE                } from './modules/annotate/parse_hmmsearch.nf'

include { GENERIC_HMM_FORMATTER                         } from './modules/annotate/generic_hmm_formatter.nf'
include { KEGG_HMM_FORMATTER                            } from './modules/annotate/kegg_hmm_formatter.nf'
include { KOFAM_HMM_FORMATTER                           } from './modules/annotate/kofam_hmm_formatter.nf'
include { DBCAN_HMM_FORMATTER                           } from './modules/annotate/dbcan_hmm_formatter.nf'
include { VOG_HMM_FORMATTER                             } from './modules/annotate/vog_hmm_formatter.nf'
include { CAMPER_HMM_FORMATTER                          } from './modules/annotate/camper_hmm_formatter.nf'
include { CANTHYD_HMM_FORMATTER                         } from './modules/annotate/canthyd_hmm_formatter.nf'
include { SULFUR_HMM_FORMATTER                          } from './modules/annotate/sulfur_hmm_formatter.nf'
include { FEGENIE_HMM_FORMATTER                          } from './modules/annotate/fegenie_hmm_formatter.nf'

include { COMBINE_ANNOTATIONS                           } from './modules/annotate/combine_annotations.nf'
include { COUNT_ANNOTATIONS                             } from './modules/annotate/count_annotations.nf'
include { MERGE_ANNOTATIONS                             } from './modules/annotate/merge_annotations.nf'

include { COMBINE_DISTILL                               } from './modules/distill/combine_distill.nf'
include { DISTILL_SUMMARY                               } from './modules/distill/distill_summary.nf'
include { DISTILL_FINAL                                 } from './modules/distill/distill_final.nf'

include { PRODUCT_HEATMAP                               } from './modules/product/product_heatmap.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Help menu and Version menu check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/* Call Help Menu */
if ((params.help || params.h) && params.call ){
    callHelpMessage()
    exit 0
}
/* Annotate Help Menu */
else if ((params.help || params.h) && params.annotate ){
    annotateHelpMessage()
    exit 0
}

/* Distill Help Menu */
else if ((params.help || params.h) && (params.distill_topic != "" || params.distill_ecosystem != "" || params.distill_custom != "") ){
    distillHelpMessage()
    exit 0
}

/* Adjectives Help Menu */
else if ((params.help || params.h) && params.adjectives ){
    adjectivesHelpMessage()
    exit 0
}

/* Options to display the current version */
else if((params.version) || (params.v)){
    version()
    exit 0
}

/* Display help message */
else if ((params.help) || (params.h)){
    helpMessage()
    exit 0
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validate Input parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def validOptions = ["--call", "--annotate", "--distill_topic", "--distill_ecosystem", "--distill_custom"]

if (params.call == 0 && params.annotate == 0 && params.annotations == "" && (params.distill_topic == "" || params.distill_ecosystem == "" || params.distill_custom == "" )) {
    error("Please provide one of the following options: ${validOptions.join(', ')}")
}

if( params.use_dbset){
    if (!['metabolism_kegg_set', 'metabolism_set', 'adjectives_kegg_set', 'adjectives_set'].contains(params.use_dbset)) {
        error("Invalid parameter '--use_dbset ${params.use_dbset}'. Valid values are 'metabolism_kegg_set', 'metabolism_set', 'adjectives_kegg_set', 'adjectives_set'.")
    }
}


if( params.annotations == "" && params.annotate == 0 && (params.distill_topic != "" || params.distill_ecosystem != "" || params.distill_custom != "" )){
    error("If you want to distill, you must provide annotations via --annotations <path/to/file>.")
}



//Add in other checks for adjectives,... etc.


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Parse DRAM2 ANNOTATE input databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Here, if --annotate is present, then parse the input databases and set variables to
// T/F for each possible database and custom database

// For example, if user uses --use_dbset metabolism, set all the individual databases
// equal to T/F (0/1)

annotate_kegg = 0
annotate_kofam = 0
annotate_pfam = 0
annotate_dbcan = 0
annotate_camper = 0
annotate_fegenie = 0
annotate_methyl = 0
annotate_canthyd = 0
annotate_heme = 0
annotate_sulfur = 0
annotate_methyl = 0
annotate_merops = 0
annotate_vogdb = 0
annotate_uniref = 0
annotate_viral = 0

if( params.use_kegg ){
    annotate_kegg = 1
}
if( params.use_kofam ){
    annotate_kofam = 1
}
if( params.use_pfam ){
    annotate_pfam = 1
}
if( params.use_dbcan ){
    annotate_dbcan = 1
}
if( params.use_camper ){
    annotate_camper = 1
}
if( params.use_fegenie){
    annotate_fegenie = 1
}
if( params.use_methyl ){
    annotate_methyl = 1
}
if( params.use_canthyd ){
    annotate_canthyd = 1
}
if( params.use_heme ){
    annotate_heme = 1
}
if( params.use_sulfur ){
    annotate_sulfur = 1
}
if( params.use_methyl ){
    annotate_methyl = 1
}
if( params.use_merops ){
    annotate_merops = 1
}
if( params.use_vog ){
    annotate_vogdb = 1
}
if( params.use_viral ){
    annotate_viral = 1
}
if( params.use_uniref ){
    annotate_uniref = 1
}

/* Metabolism Database sets */
if( params.use_dbset == "metabolism_kegg_set" ){
    annotate_kegg = 1
    annotate_dbcan = 1
    annotate_merops = 1
    annotate_pfam = 1
    annotate_heme = 1
}
if( params.use_dbset == "metabolism_set" ){
    annotate_kofam = 1
    annotate_dbcan = 1
    annotate_merops = 1
    annotate_pfam = 1
    annotate_heme = 1
}
if( params.use_dbset == "adjectives_kegg_set" ){
    annotate_kegg = 1
    annotate_dbcan = 1
    annotate_merops = 1
    annotate_pfam = 1
    annotate_heme = 1
    annotate_sulfur = 1
    annotate_camper = 1
    annotate_methyl = 1   
    annotate_fegenie = 1
}
if( params.use_dbset == "adjectives_set" ){
    annotate_kofam = 1
    annotate_dbcan = 1
    annotate_merops = 1
    annotate_pfam = 1
    annotate_heme = 1
    annotate_sulfur = 1
    annotate_camper = 1
    annotate_methyl = 1   
    annotate_fegenie = 1   
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Create channel for reference databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//This is where we will check, for example, if --annotate is provided and if it is, create 
//  channels for the various annotation databases
    def annotate_list = "" 

if( params.annotate ){
    //This is just temporary - want these in the containers eventually
    ch_combine_annot_script = file(params.combine_annotations_script)
    ch_count_annots_script = file(params.count_annots_script)
    ch_parse_hmmsearch = file(params.parse_hmmsearch_script)
    ch_mmseqs_script = file(params.mmseqs_add_descriptions_script)

    ch_kegg_formatter = file(params.kegg_formatter_script)
    ch_kofam_formatter = file(params.kofam_hmm_formatter_script)
    ch_dbcan_formatter = file(params.dbcan_hmm_formatter_script)
    ch_vog_formatter = file(params.vog_hmm_formatter_script)
    ch_camper_formatter = file(params.camper_hmm_formatter_script)
    ch_canthyd_formatter = file(params.canthyd_hmm_formatter_script)
    ch_sulfur_formatter = file(params.sulfur_hmm_formatter_script)
    ch_fegenie_formatter = file(params.sulfur_hmm_formatter_script)

    ch_kofam_list = file(params.kofam_list)
    ch_canthyd_list = file(params.cant_hyd_hmm_list)
    ch_dbcan_fam = file(params.dbcan_fam_activities)
    ch_dbcan_subfam = file(params.dbcan_subfam_activities)
    ch_vog_list = file(params.vog_list)
    ch_camper_hmm_list = file(params.camper_hmm_list)
    ch_canthyd_hmm_list = file(params.cant_hyd_hmm_list)

    ch_dummy_sheet = file(params.distill_dummy_sheet)

    ch_sql_parser = file(params.sql_parser_script)

    ch_sql_descriptions_db = file(params.sql_descriptions_db)

    index_mmseqs = "0"

    if (annotate_kegg == 1) {
        ch_kegg_db = file(params.kegg_db).exists() ? file(params.kegg_db) : error("Error: If using --annotate, you must supply prebuilt databases. KEGG database file not found at ${params.kegg_db}")
        index_mmseqs = "1"
        annotate_list += "KEGG "
    }

    if (annotate_kofam == 1) {
        ch_kofam_db = file(params.kofam_db).exists() ? file(params.kofam_db) : error("Error: If using --annotate, you must supply prebuilt databases. KOFAM database file not found at ${params.kofam_db}")
        annotate_list += "Kofam "
    }

    if (annotate_dbcan == 1) {
        ch_dbcan_db = file(params.dbcan_db).exists() ? file(params.dbcan_db) : error("Error: If using --annotate, you must supply prebuilt databases. DBCAN database file not found at ${params.dbcan_db}")
        annotate_list += "dbCAN "
    }

    if (annotate_camper == 1) {
        ch_camper_hmm_db = file(params.camper_hmm_db).exists() ? file(params.camper_hmm_db) : error("Error: If using --annotate, you must supply prebuilt databases. CAMPER HMM database file not found at ${params.camper_hmm_db}")
        ch_camper_mmseqs_db = file(params.camper_mmseqs_db).exists() ? file(params.camper_mmseqs_db) : error("Error: If using --annotate, you must supply prebuilt databases. CAMPER MMseqs2 database file not found at ${params.camper_mmseqs_db}")
        index_mmseqs = "1"
        annotate_list += "CAMPER "
        ch_camper_mmseqs_list = file(params.camper_mmseqs_list)
    }

    if (annotate_merops == 1) {
        ch_merops_db = file(params.merops_db).exists() ? file(params.merops_db) : error("Error: If using --annotate, you must supply prebuilt databases. MEROPS database file not found at ${params.merops_db}")
        index_mmseqs = "1"
        annotate_list += "MEROPS "
    }

    if (annotate_pfam == 1) {
        ch_pfam_mmseqs_db = file(params.pfam_mmseq_db).exists() ? file(params.pfam_mmseq_db) : error("Error: If using --annotate, you must supply prebuilt databases. PFAM database file not found at ${params.pfam_mmseq_db}")
        index_mmseqs = "1"
        annotate_list += "Pfam "
    }

    if (annotate_heme == 1) {
        ch_heme_db = file(params.heme_db).exists() ? file(params.heme_db) : error("Error: If using --annotate, you must supply prebuilt databases. HEME database file not found at ${params.heme_db}")
        annotate_list += "hene "
    }

    if (annotate_sulfur == 1) {
        ch_sulfur_db = file(params.sulfur_db).exists() ? file(params.sulfur_db) : error("Error: If using --annotate, you must supply prebuilt databases. SULURR database file not found at ${params.sulfur_db}")
        annotate_list += "sulfur "
    }

    if (annotate_uniref == 1) {
        ch_uniref_db = file(params.uniref_db).exists() ? file(params.uniref_db) : error("Error: If using --annotate, you must supply prebuilt databases. UNIREF database file not found at ${params.uniref_db}")
        index_mmseqs = "1"
        annotate_list += "UniRef "
    }

    if (annotate_methyl == 1) {
        ch_methyl_db = file(params.methyl_db).exists() ? file(params.methyl_db) : error("Error: If using --annotate, you must supply prebuilt databases. METHYL database file not found at ${params.methyl_db}")
        index_mmseqs = "1"
        annotate_list += "methyl "
    }

    if (annotate_fegenie == 1) {
        ch_fegenie_db = file(params.fegenie_db).exists() ? file(params.fegenie_db) : error("Error: If using --annotate, you must supply prebuilt databases. FEGENIE database file not found at ${params.fegenie_db}")
        annotate_list += "FeGenie "
    }

    if (annotate_canthyd == 1) {
        ch_canthyd_hmm_db = file(params.canthyd_hmm_db).exists() ? file(params.canthyd_hmm_db) : error("Error: If using --annotate, you must supply prebuilt databases. CANT_HYD HMM database file not found at ${params.canthyd_hmm_db}")
        ch_canthyd_mmseqs_db = file(params.canthyd_mmseqs_db).exists() ? file(params.canthyd_mmseqs_db) : error("Error: If using --annotate, you must supply prebuilt databases. CANT_HYD MMseqs database file not found at ${params.canthyd_mmseqs_db}")
        index_mmseqs = "1"
        annotate_list += "CANT-HYD "
        ch_canthyd_mmseqs_list = file(params.canthyd_mmseqs_list)
    }

    if (annotate_vogdb == 1) {
        ch_vogdb_db = file(params.vog_db).exists() ? file(params.vog_db) : error("Error: If using --annotate, you must supply prebuilt databases. VOG database file not found at ${params.vog_db}")
        annotate_list += "VOGDB "
    }

    if (annotate_viral == 1) {
        ch_viral_db = file(params.viral_db).exists() ? file(params.viral_db) : error("Error: If using --annotate, you must supply prebuilt databases. viral database file not found at ${params.viral_db}")
        index_mmseqs = "1"
        annotate_list += "viral "
    }
    /* Custom user databases */
    //Not sure about this yet.
    // Consider: hmm format, mmseqs2 format....

    // Continue for the rest of the databases

    /* Check for input Bin Quality file */
    if (params.bin_quality != "") {
        ch_bin_quality = file(params.bin_quality).exists() ? file(params.bin_quality) : error("Error: If using --bin_quality, you must supply a formatted input file. Bin quality file not found at ${params.bin_quality}")
    } else {
        ch_bin_quality = []
    }

    /* Check for input Taxa file */
    if (params.taxa != "") {
        ch_taxa = file(params.taxa).exists() ? file(params.taxa) : error("Error: If using --taxa, you must supply a formatted input file. Taxonomy file not found at ${params.taxa}")
    } else {
        ch_taxa = []
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Create channel for INGEST inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def distill_flag = "off"

// Here we will check for the various kinds of inputs
// The structure of this will change
// For example, we will be able to ingest only metaT data and call genes on that
// For example, users may use a csv or put a path to a directory
if( params.call ){
    // If calling genes, then create a channel called ch_input_fastas.
    if ( params.input_fasta ) {
        ch_input_fastas = Channel
            .fromPath(params.input_fasta + params.fasta_fmt, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find any fasta files matching: ${params.input_fasta}\nNB: Path needs to follow pattern: path/to/directory/" }
    } else {
        ch_input_fastas = Channel
            .fromPath(params.fastas, checkIfExists: true)
            .ifEmpty { exit 1, "If you specify --input_fasta you must provide a path/to/directory containing fasta files or they must be in: ./raw_data/*.fa" }
    }    

    // Validate prodigal options for Call
    if (!['single', 'meta'].contains(params.prodigal_mode)) {
        error("Invalid parameter '--prodigal_mode ${params.prodigal_mode}'. Valid values are 'single' or 'meta'.")
    }
    if (!['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25'].contains(params.prodigal_trans_table)) {
        error("Invalid parameter '--prodigal_trans_table ${params.binning_map_mode}'. Valid values are '1','2',...,'25'.")
    }

    // Convert ch_input_fastas into a tuple: [samplename, file]
    ch_input_fastas.map {
        sampleName = it.getName().replaceAll(/\.[^.]+$/, '').replaceAll(/\./, '-')
        tuple(sampleName, it)
    }.set{ch_input_fastas}
}

if( params.annotate ){
    /* If the user did not specify --call, then set input called genes and proteins */
    if( params.call == 0 ){
        // Set ch_input_genes
        ch_called_genes = Channel
            .fromPath(params.input_genes + params.genes_fmt, checkIfExists: true)
            .ifEmpty { exit 1, "If you specify --annotate without --call, you must provide a fasta file of called genes. Cannot find any called gene fasta files matching: ${params.input_genes}\nNB: Path needs to follow pattern: path/to/directory/" }
            .map {
                sampleName = it.getName().replaceAll(/\.[^.]+$/, '').replaceAll(/\./, '-')
                tuple(sampleName, it)
            }

        // Set ch_input_proteins
        ch_called_proteins = Channel
            .fromPath(params.input_proteins + params.proteins_fmt, checkIfExists: true)
            .ifEmpty { exit 1, "If you specify --annotate without --call, you must provide a fasta file of called proteins. Cannot find any called gene fasta files matching: ${params.input_proteins}\nNB: Path needs to follow pattern: path/to/directory/" }
            .map {
                sampleName = it.getName().replaceAll(/\.[^.]+$/, '').replaceAll(/\./, '-')
                tuple(sampleName, it)
            }
    }    
    
    /* Check for input Bin Quality file */
    if (params.bin_quality != "") {
        ch_bin_quality = file(params.bin_quality).exists() ? file(params.bin_quality) : error("Error: If using --bin_quality, you must supply a formatted input file. Bin quality file not found at ${params.bin_quality}")
    } else {
        ch_bin_quality = []
    }

    /* Check for input Taxa file */
    if (params.taxa != "") {
        ch_taxa = file(params.taxa).exists() ? file(params.taxa) : error("Error: If using --taxa, you must supply a formatted input file. Taxonomy file not found at ${params.taxa}")
    } else {
        ch_taxa = []
    } 

    // Ensure an add_annotations channel is generated if the user specifies --add_annotations
    if( params.add_annotations != ""){
        ch_add_annots = file(params.add_annotations).exists() ? file(params.add_annotations) : error("Error: If using --add_annotations, you must supply a DRAM-formatted annotations file. Taxonomy file not found at ${params.add_annotations}")
    }


}

if (params.distill_topic != "" || params.distill_ecosystem != "" || params.distill_custom != "") { 
    distill_flag = "on"
    // Set channels for supporting python scripts - will be moved to container eventually
    ch_distill_summary_script = file(params.distill_summary_script)
    ch_distill_final_script = file(params.distill_final_script)  

    // Ensure rRNA and tRNA channels are populated if the user is not calling genes
    if( params.call == 0 ){
        //set ch_rrna_sheet = RRNA_COLLECT.out.rrna_collected_out
        //    ch_rrna_combined = RRNA_COLLECT.out.rrna_combined_out
        //    ch_trna_sheet = TRNA_COLLECT.out.trna_collected_out
        if( params.rrnas != "" ){
            Channel.fromPath("${params.rrnas}/*.tsv", checkIfExists: true)
                .ifEmpty { exit 1, "If you specify --distill_<topic|ecosystem|custom> without --call, you must provide individual rRNA files generated with barrnap. Cannot find any files at: ${params.rrnas}\nNB: Path needs to follow pattern: path/to/directory" }
                .collect()
                .set { ch_collected_rRNAs }
        }
        if( params.trnas != "" ){
            Channel.fromPath("${params.trnas}/*.tsv", checkIfExists: true)
                .ifEmpty { exit 1, "If you specify --distill_<topic|ecosystem|custom> without --call, you must provide individual rRNA files generated with tRNAscan-SE. Cannot find any files at: ${params.trnas}\nNB: Path needs to follow pattern: path/to/directory" }
                .collect()
                .set { ch_collected_tRNAs }
        }
    }

    // Ensure annotations, taxonomy and bin quality channels are set.
    if( params.annotate == 0 ){
        // Set channels for of supporting python scripts - will be moved to container eventually
        ch_count_annots_script = file(params.count_annots_script)

        ch_combined_annotations = Channel
            .fromPath(params.annotations, checkIfExists: true)
            .ifEmpty { exit 1, "If you specify --distill_<topic|ecosystem|custom> without --annotate, you must provide an annotations TSV file (--annotations <path>) with approprite formatting. Cannot find any called gene fasta files matching: ${params.annotations}\nNB: Path needs to follow pattern: path/to/directory/" }
        
        /* Check for input Bin Quality file */
        if (params.bin_quality != "") {
            ch_bin_quality = file(params.bin_quality).exists() ? file(params.bin_quality) : error("Error: If using --bin_quality, you must supply a formatted input file. Bin quality file not found at ${params.bin_quality}")
        } else {
            ch_bin_quality = params.distill_dummy_sheet
        }

        /* Check for input Taxa file */
        if (params.taxa != "") {
            ch_taxa = file(params.taxa).exists() ? file(params.taxa) : error("Error: If using --taxa, you must supply a formatted input file. Taxonomy file not found at ${params.taxa}")
        } else {
            ch_taxa = params.distill_dummy_sheet
        }

        // Ensure an add_annotations channel is generated if the user specifies --add_annotations
        if( params.add_annotations != ""){
            ch_add_annots = file(params.add_annotations).exists() ? file(params.add_annotations) : error("Error: If using --add_annotations, you must supply a DRAM-formatted annotations file. Taxonomy file not found at ${params.add_annotations}")
        }

    }



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Create channels for optional topic and/or ecosystem and/or custom sheets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/* Create the default distill topic and ecosystem channels */
default_channel = Channel.fromPath(params.distill_dummy_sheet)
def distill_topic_list = "" 
def distill_ecosystem_list = ""
def distill_custom_list = ""

distill_default = "0"
distill_carbon = "0"
distill_energy = "0"
distill_misc = "0"
distill_nitrogen = "0"
distill_transport = "0"

if (params.distill_topic != "" || params.distill_ecosystem != "" || params.distill_custom != "") {    
    if (params.distill_topic != "") {
        def validTopics = ['default', 'carbon', 'energy', 'misc', 'nitrogen', 'transport']
        def topics = params.distill_topic.split()
        
        // If 'default' is one of the topics, ignore other topics and set everything to "1"
        if (topics.contains("default")) {
            distill_default = "1"
            distill_carbon = "1"
            distill_energy = "1"
            distill_misc = "1"
            distill_nitrogen = "1"
            distill_transport = "1"
            distill_topic_list = "default (carbon, energy, misc, nitrogen, transport)"
        } else {
            // Process other topics only if 'default' is not selected
            topics.each { topic ->
                if (!validTopics.contains(topic)) {
                    error("Invalid distill topic: $topic. Valid values are ${validTopics.join(', ')}")
                }

                switch (topic) {
                    case "carbon":
                        distill_carbon = "1"
                        distill_topic_list += "carbon "
                        break
                    case "energy":
                        distill_energy = "1"
                        distill_topic_list += "energy "
                        break
                    case "misc":
                        distill_misc = "1"
                        distill_topic_list += "misc "
                        break
                    case "nitrogen":
                        distill_nitrogen = "1"
                        distill_topic_list += "nitrogen "
                        break
                    case "transport":
                        distill_transport = "1"
                        distill_topic_list += "transport "
                        break
                }
            }
        }
    }

    if (distill_carbon == "1") {
        ch_distill_carbon = file(params.distill_carbon_sheet).exists() ? file(params.distill_carbon_sheet) : error("Error: If using --distill_topic carbon (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

    } else{
        ch_distill_carbon = default_channel
    }
    if (distill_energy == "1") {
        ch_distill_energy = file(params.distill_energy_sheet).exists() ? file(params.distill_energy_sheet) : error("Error: If using --distill_topic energy (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

    } else{
        ch_distill_energy = default_channel
    }

    if (distill_misc == "1") {
        ch_distill_misc = file(params.distill_misc_sheet).exists() ? file(params.distill_misc_sheet) : error("Error: If using --distill_topic misc (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

    } else{
        ch_distill_misc = default_channel
    }

    if (distill_nitrogen == "1") {
        ch_distill_nitrogen = file(params.distill_nitrogen_sheet).exists() ? file(params.distill_nitrogen_sheet) : error("Error: If using --distill_topic nitrogen (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

    } else{
        ch_distill_nitrogen = default_channel
    }

    if (distill_transport == "1") {
        ch_distill_transport = file(params.distill_transport_sheet).exists() ? file(params.distill_transport_sheet) : error("Error: If using --distill_topic transport (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

    } else{
        ch_distill_transport = default_channel
    }

    distill_eng_sys = "0"
    distill_ag = "0"   
    if (params.distill_ecosystem != "") {
        def distillEcosystemList = params.distill_ecosystem.split()

        // Create a list to store the generated channels
        def ecoSysChannels = []

        distillEcosystemList.each { ecosysItem ->
            if (!['eng_sys', 'ag'].contains(ecosysItem)) {
                error("Invalid distill ecosystem: $ecosysItem. Valid values are eng_sys, ag")
            }

            switch (ecosysItem) {
                case "ag":
                    distill_ag = "1"
                    distill_ecosystem_list = distill_ecosystem_list + "ag "
                    break
                case "eng_sys":
                    distill_eng_sys = "1"
                    distill_ecosystem_list = distill_ecosystem_list + "eng_sys "
                    break
            }
        }
    }

    if (distill_eng_sys == "1") {
        ch_distill_eng_sys = file(params.distill_eng_sys_sheet).exists() ? file(params.distill_eng_sys_sheet) : error("Error: If using --distill_ecosystem eng_sys, you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

    } else{
        ch_distill_eng_sys = default_channel
    }

    if (distill_ag == "1") {
        ch_distill_ag = file(params.distill_ag_sheet).exists() ? file(params.distill_ag_sheet) : error("Error: If using --distill_ecosystem ag, you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

    } else{
        ch_distill_ag = default_channel
    }

    if (params.distill_custom != "") {
        ch_distill_custom = file(params.distill_custom).exists() ? file(params.distill_custom) : error("Error: If using --distill_custom <path/to/TSV>, you must have the preformatted custom distill sheet in the provided file: ${params.distill_custom}.")
        distill_custom_list = "params.distill_custom"
    }else{
        ch_distill_custom = default_channel
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Print run info to command line
        Various run info for various DRAM2 options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// This is just a catch-all for now - NEED to generate others for various options
if( params.call && params.annotate && (params.distill_ecosystem !="" || params.distill_custom !="" || params.distill_topic !="" )){
    log.info """
            DRAM2 Nextflow
            ===================================
            fastas       : ${params.input_fasta}
            outdir       : ${params.outdir}
            threads      : ${params.threads}
            rename       : ${params.rename ? 'true' : 'false'}
            call genes   : ${params.call ? 'true' : 'false'}
            annotate     : ${params.annotate ? 'true' : 'false'}
            databases    : ${annotate_list}
            distill      : ${distill_flag}
              topic      : ${distill_topic_list}
              ecosystem  : ${distill_ecosystem_list}
              custom     : ${params.distill_custom}

            """
            .stripIndent()
}else if( !params.call && params.annotate && (params.distill_ecosystem !="" || params.distill_custom !="" || params.distill_topic !="" )){
    log.info """
            DRAM2 Nextflow
            ===================================
            fastas       : ${params.input_fasta}
            outdir       : ${params.outdir}
            threads      : ${params.threads}
            annotate     : ${params.annotate ? 'true' : 'false'}
            databases    : ${annotate_list}
            distill      : ${distill_flag}
              topic      : ${distill_topic_list}
              ecosystem  : ${distill_ecosystem_list}
              custom     : ${params.distill_custom}

            """
            .stripIndent()
}else if( params.call == 0 && params.annotate == 0 && (params.distill_ecosystem !="" || params.distill_custom !="" || params.distill_topic !="" )){
    log.info """
            DRAM2 Nextflow
            ===================================
            annotations  : ${params.annotations}
            outdir       : ${params.outdir}
            threads      : ${params.threads}
            tRNA         : ${params.trnas}
            rRNA         : ${params.rrnas}
            databases    : 
            distill      : ${distill_flag}
              topic      : ${distill_topic_list}
              ecosystem  : ${distill_ecosystem_list}
              custom     : ${params.distill_custom}

            """
            .stripIndent()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Main workflow for pipeline

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    /* Rename fasta headers
        Process 1-by-1 
    */
    if( params.call ){
        if( params.rename ) {
            RENAME_FASTA( ch_input_fastas )
            ch_fasta = RENAME_FASTA.out.renamed_fasta
        } 
        else {
            ch_fasta = ch_input_fastas
        }

        CALL_GENES ( ch_fasta )
        called_genes = CALL_GENES.out.prodigal_fna
        ch_called_proteins = CALL_GENES.out.prodigal_faa

        // Not sure if we want these default or optional yet
        /* Run tRNAscan-SE on each fasta to identify tRNAs */
        TRNA_SCAN( ch_fasta )
        ch_trna_scan = TRNA_SCAN.out.trna_scan_out
        // Collect all sample formatted tRNA files
        Channel.empty()
            .mix( ch_trna_scan )
            .collect()
            .set { ch_collected_tRNAs }
        /* Run TRNA_COLLECT to generate a combined TSV for all fastas */
        TRNA_COLLECT( ch_collected_tRNAs )
        ch_trna_sheet = TRNA_COLLECT.out.trna_collected_out

        /* Run barrnap on each fasta to identify rRNAs */
        RRNA_SCAN( ch_fasta )
        ch_rrna_scan = RRNA_SCAN.out.rrna_scan_out
        Channel.empty()
            .mix( ch_rrna_scan )
            .collect()
            .set { ch_collected_rRNAs }
        /* Run RRNA_COLLECT to generate a combined TSV for all fastas */
        RRNA_COLLECT( ch_collected_rRNAs )
        ch_rrna_sheet = RRNA_COLLECT.out.rrna_collected_out
        ch_rrna_combined = RRNA_COLLECT.out.rrna_combined_out


    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Annotation
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if( params.annotate ){

        def formattedOutputChannels = channel.of()

        // Here we will create mmseqs2 index files for each of the inputs if we are going to do a mmseqs2 database
        if( index_mmseqs == "1" ){
            MMSEQS_INDEX( ch_called_proteins )
            ch_mmseqs_query = MMSEQS_INDEX.out.mmseqs_index_out
        }

        // Annotate according to the user-specified databases 
        if( annotate_kegg == 1 ){
            MMSEQS_SEARCH_KEGG( ch_mmseqs_query, ch_kegg_db, params.bit_score_threshold, ch_dummy_sheet, params.kegg_name, ch_mmseqs_script )
            ch_kegg_unformatted = MMSEQS_SEARCH_KEGG.out.mmseqs_search_formatted_out

            SQL_KEGG(ch_kegg_unformatted, params.kegg_name, ch_sql_descriptions_db, ch_sql_parser)
            ch_kegg_formatted = SQL_KEGG.out.sql_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_kegg_formatted)
        }

        if( annotate_kofam == 1 ){
            HMM_SEARCH_KOFAM ( ch_called_proteins,  params.kofam_e_value, ch_kofam_db )
            ch_kofam_hmms = HMM_SEARCH_KOFAM.out.hmm_search_out

            PARSE_HMM_KOFAM ( ch_kofam_hmms, ch_parse_hmmsearch )
            ch_kofam_parsed = PARSE_HMM_KOFAM.out.parsed_hmm

            KOFAM_HMM_FORMATTER ( ch_kofam_parsed, params.kofam_top_hit, ch_kofam_list, ch_kofam_formatter )
            ch_kofam_formatted = KOFAM_HMM_FORMATTER.out.kofam_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_kofam_formatted)
        }
        //NOT DONE
        if( annotate_pfam == 1 ){
            MMSEQS_SEARCH_PFAM( ch_mmseqs_query, ch_pfam_mmseqs_db, params.bit_score_threshold, ch_dummy_sheet, params.pfam_name, ch_mmseqs_script )
            ch_uniref_formatted = MMSEQS_SEARCH_PFAM.out.mmseqs_search_formatted_out

            formattedOutputChannels = formattedOutputChannels.mix(ch_uniref_formatted)
        }

        if( annotate_dbcan == 1 ){
            
            HMM_SEARCH_DBCAN ( ch_called_proteins, params.dbcan_e_value , ch_dbcan_db)
            ch_dbcan_hmms = HMM_SEARCH_DBCAN.out.hmm_search_out

            PARSE_HMM_DBCAN ( ch_dbcan_hmms, ch_parse_hmmsearch )
            ch_dbcan_parsed = PARSE_HMM_DBCAN.out.parsed_hmm

            DBCAN_HMM_FORMATTER ( ch_dbcan_parsed, params.dbcan_top_hit, params.dbcan_name, ch_dbcan_formatter, ch_sql_parser, ch_sql_descriptions_db )
            ch_dbcan_formatted = DBCAN_HMM_FORMATTER.out.dbcan_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_dbcan_formatted)
        }

        if (annotate_camper == 1){
            // HMM
            HMM_SEARCH_CAMPER ( ch_called_proteins, params.camper_e_value , ch_camper_hmm_db)
            ch_camper_hmms = HMM_SEARCH_CAMPER.out.hmm_search_out

            PARSE_HMM_CAMPER ( ch_camper_hmms, ch_parse_hmmsearch )
            ch_camper_parsed = PARSE_HMM_CAMPER.out.parsed_hmm

            CAMPER_HMM_FORMATTER ( ch_camper_parsed, params.camper_top_hit, ch_camper_hmm_list, ch_camper_formatter )
            ch_camper_hmm_formatted = CAMPER_HMM_FORMATTER.out.camper_formatted_hits
            
            formattedOutputChannels = formattedOutputChannels.mix(ch_camper_hmm_formatted)

            // MMseqs
            MMSEQS_SEARCH_CAMPER( ch_mmseqs_query, ch_camper_mmseqs_db, params.bit_score_threshold, ch_camper_mmseqs_list, params.camper_name, ch_mmseqs_script )
            ch_camper_mmseqs_formatted = MMSEQS_SEARCH_CAMPER.out.mmseqs_search_formatted_out

            formattedOutputChannels = formattedOutputChannels.mix(ch_camper_mmseqs_formatted)
        }
        // NOT DONE - HMM
        if (annotate_fegenie == 1){
            HMM_SEARCH_FEGENIE ( ch_called_proteins,  params.fegenie_e_value, ch_fegenie_db )
            ch_fegenie_hmms = HMM_SEARCH_FEGENIE.out.hmm_search_out

            PARSE_HMM_FEGENIE ( ch_fegenie_hmms, ch_parse_hmmsearch )
            ch_fegenie_parsed = PARSE_HMM_FEGENIE.out.parsed_hmm

            FEGENIE_HMM_FORMATTER ( ch_fegenie_parsed, ch_fegenie_formatter )
            ch_fegenie_formatted = FEGENIE_HMM_FORMATTER.out.fegenie_formatted_hits
            formattedOutputChannels = formattedOutputChannels.mix(ch_fegenie_formatted)
        }

        if (annotate_methyl == 1){
            MMSEQS_SEARCH_METHYL( ch_mmseqs_query, ch_methyl_db, params.bit_score_threshold, ch_dummy_sheet, params.methyl_name, ch_mmseqs_script )
            ch_methyl_mmseqs_formatted = MMSEQS_SEARCH_METHYL.out.mmseqs_search_formatted_out

            formattedOutputChannels = formattedOutputChannels.mix(ch_methyl_mmseqs_formatted)
        }

        if (annotate_canthyd == 1){
            // MMseqs
            MMSEQS_SEARCH_CANTHYD( ch_mmseqs_query, ch_canthyd_mmseqs_db, params.bit_score_threshold, ch_canthyd_mmseqs_list, params.canthyd_name, ch_mmseqs_script )
            ch_canthyd_mmseqs_formatted = MMSEQS_SEARCH_CANTHYD.out.mmseqs_search_formatted_out

            formattedOutputChannels = formattedOutputChannels.mix(ch_canthyd_mmseqs_formatted)

            //HMM
            HMM_SEARCH_CANTHYD ( ch_called_proteins, params.canthyd_e_value , ch_canthyd_hmm_db)
            ch_canthyd_hmms = HMM_SEARCH_CANTHYD.out.hmm_search_out

            PARSE_HMM_CANTHYD ( ch_canthyd_hmms, ch_parse_hmmsearch )
            ch_canthyd_parsed = PARSE_HMM_CANTHYD.out.parsed_hmm

            CANTHYD_HMM_FORMATTER ( ch_canthyd_parsed, params.canthyd_top_hit, ch_canthyd_hmm_list, ch_canthyd_formatter )
            ch_canthyd_hmm_formatted = CANTHYD_HMM_FORMATTER.out.canthyd_formatted_hits
            
            formattedOutputChannels = formattedOutputChannels.mix(ch_canthyd_hmm_formatted)

        }
        // Not done - not sure what this is?
        if (annotate_heme == 1){
            formattedOutputChannels =  formattedOutputChannels.mix(ch_heme_formatted)
        }

        if (annotate_sulfur == 1){
            HMM_SEARCH_SULFUR ( ch_called_proteins,  params.sulfur_e_value, ch_sulfur_db )
            ch_sulfur_hmms = HMM_SEARCH_SULFUR.out.hmm_search_out

            PARSE_HMM_SULFUR ( ch_sulfur_hmms, ch_parse_hmmsearch )
            ch_sulfur_parsed = PARSE_HMM_SULFUR.out.parsed_hmm

            SULFUR_HMM_FORMATTER ( ch_sulfur_parsed, ch_sulfur_formatter )
            ch_sulfur_formatted = SULFUR_HMM_FORMATTER.out.sulfur_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_sulfur_formatted)
        }

        if (annotate_merops == 1){
            MMSEQS_SEARCH_MEROPS( ch_mmseqs_query, ch_merops_db, params.bit_score_threshold, ch_dummy_sheet, params.merops_name, ch_mmseqs_script )
            ch_merops_unformatted = MMSEQS_SEARCH_MEROPS.out.mmseqs_search_formatted_out

            SQL_MEROPS(ch_merops_unformatted, params.merops_name, ch_sql_descriptions_db, ch_sql_parser)
            ch_merops_formatted = SQL_MEROPS.out.sql_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_merops_formatted)
        }

        if (annotate_uniref == 1){
            MMSEQS_SEARCH_UNIREF( ch_mmseqs_query, ch_uniref_db, params.bit_score_threshold, ch_dummy_sheet, params.uniref_name, ch_mmseqs_script )
            ch_uniref_unformatted = MMSEQS_SEARCH_UNIREF.out.mmseqs_search_formatted_out

            SQL_UNIREF(ch_uniref_unformatted, params.uniref_name, ch_sql_descriptions_db, ch_sql_parser)
            ch_uniref_formatted = SQL_UNIREF.out.sql_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_uniref_formatted)
        }

        if (annotate_vogdb == 1){
            HMM_SEARCH_VOG ( ch_called_proteins, params.vog_e_value , ch_vogdb_db )
            ch_vog_hmms = HMM_SEARCH_VOG.out.hmm_search_out

            PARSE_HMM_VOG ( ch_vog_hmms, ch_parse_hmmsearch )
            ch_vog_parsed = PARSE_HMM_VOG.out.parsed_hmm

            VOG_HMM_FORMATTER ( ch_vog_parsed, params.vog_top_hit, params.vogdb_name, ch_vog_formatter, ch_sql_parser, ch_sql_descriptions_db )
            ch_vog_formatted = VOG_HMM_FORMATTER.out.vog_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_vog_formatted)
        }

        if (annotate_viral == 1){
            MMSEQS_SEARCH_VIRAL( ch_mmseqs_query, ch_viral_db, params.bit_score_threshold, ch_dummy_sheet, params.viral_name, ch_mmseqs_script )
            ch_viral_unformatted = MMSEQS_SEARCH_VIRAL.out.mmseqs_search_formatted_out

            SQL_VIRAL(ch_viral_unformatted, params.viral_name, ch_sql_descriptions_db, ch_sql_parser)
            ch_viral_formatted = SQL_UNIREF.out.sql_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_viral_formatted)
        }

        Channel.empty()
            .mix( formattedOutputChannels )
            .collect()
            .set { collected_formatted_hits }

        // COMBINE_ANNOTATIONS collects all annotations files across ALL databases 
        COMBINE_ANNOTATIONS( collected_formatted_hits, ch_combine_annot_script )
        ch_combined_annotations = COMBINE_ANNOTATIONS.out.combined_annotations_out

        // Add Bin Quality to annotations 
        if( params.bin_quality != "" ){
            ADD_BIN_QUALITY( ch_combined_annotations, ch_bin_quality )
            ch_updated_annots = ADD_BIN_QUALITY.out.annots_bin_quality_out
        }
        else{
            ch_updated_annots = ch_combined_annotations
        }

        // Add Taxonomy to annotations 
        if( params.taxa != "" ){
            ADD_TAXA( ch_updated_annots, ch_taxa )
            ch_updated_taxa_annots = ADD_TAXA.out.annots_taxa_out
        }
        else{
            ch_updated_taxa_annots = ch_combined_annotations
        }

        // Check for additional user-provided annotations 
        if( params.add_annotations != "" ){
            MERGE_ANNOTATIONS( ch_updated_taxa_annots, ch_add_annots )
            ch_final_annots = MERGE_ANNOTATIONS.out.merged_annots_out
            
            COUNT_ANNOTATIONS ( ch_final_annots, ch_count_annots_script )
            ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts
        }
        else{
            ch_final_annots = ch_combined_annotations
            COUNT_ANNOTATIONS ( ch_final_annots, ch_count_annots_script )
            ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts
        }

        

    }
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Distill
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */   
    if( params.distill_topic != "" || params.distill_ecosystem != "" || params.distill_custom != "" )
    {
        // If the user did not call genes, collect tRNA and rRNA individual files provided by --trnas and --rrnas
        if( params.call == 0 ){\
            if( params.trnas != "" ){
                TRNA_COLLECT( ch_collected_tRNAs )
                ch_trna_sheet = TRNA_COLLECT.out.trna_collected_out
            }else{
                ch_trna_sheet = params.distill_dummy_sheet
            }
            if( params.rrnas != "" ){
            RRNA_COLLECT( ch_collected_rRNAs )
                ch_rrna_sheet = RRNA_COLLECT.out.rrna_collected_out
                ch_rrna_combined = RRNA_COLLECT.out.rrna_combined_out
            }else{
                ch_rrna_sheet = params.distill_dummy_sheet
                ch_rrna_combined = params.distill_dummy_sheet
            }
        }       
        // If the user did not annotate and provided taxa and/or bin quality, add it to annotations.
        if( params.annotate == 0 ){
            /* Add Bin Quality to annotations */
            if( params.bin_quality != "" ){
                ADD_BIN_QUALITY( ch_combined_annotations, ch_bin_quality )
                ch_updated_annots = ADD_BIN_QUALITY.out.annots_bin_quality_out
            }
            else{
                ch_updated_annots = ch_combined_annotations
            }
            /* Add Taxonomy to annotations */
            if( params.taxa != "" ){
                ADD_TAXA( ch_updated_annots, ch_taxa )
                ch_updated_taxa_annots = ADD_TAXA.out.annots_taxa_out
            }
            else{
                ch_updated_taxa_annots = ch_combined_annotations
            }
            if( params.add_annotations != "" ){
                MERGE_ANNOTATIONS( ch_updated_taxa_annots, ch_add_annots )
                ch_final_annots = MERGE_ANNOTATIONS.out.merged_annots_out
                
                COUNT_ANNOTATIONS ( ch_final_annots, ch_count_annots_script )
                ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts
            }
            else{
                ch_final_annots = ch_combined_annotations
                COUNT_ANNOTATIONS ( ch_final_annots, ch_count_annots_script )
            }
        } 
        
    
        
        /* Combine the individual user-specified distill sheets into a single channel */
        COMBINE_DISTILL(ch_distill_carbon, ch_distill_energy, ch_distill_misc, ch_distill_nitrogen, ch_distill_transport, ch_distill_ag, ch_distill_eng_sys, ch_distill_custom )
        ch_combined_distill_sheets = COMBINE_DISTILL.out.ch_combined_distill_sheets

        /* Generate a single distillate sheet which will then be separated by DISTILL_FINAL */
        DISTILL_SUMMARY( ch_final_annots, ch_combined_distill_sheets, ch_annotation_counts, ch_distill_summary_script )
        ch_simple_matab_summ = DISTILL_SUMMARY.out.ch_genome_sum_simple

        /* Separate the distill summary into separate sheets - add on genome_stats sheet, rRNA sheet and tRNA sheet */
        DISTILL_FINAL( ch_simple_matab_summ, ch_distill_final_script, ch_rrna_sheet, ch_rrna_combined, ch_trna_sheet, ch_final_annots )

        
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Product
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */   
    /*
    if( params.product ){
        PRODUCT_HEATMAP( ch_annotation_counts )

    }
    */

}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Version menu
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def version() {
    log.info"""

    DRAM2 

    Software versions used:
    BBTools    v39.01
    Bowtie2    v2.5.1

    """.stripIndent()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Help menus
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def helpMessage() {
    log.info """
    DRAM2 Nextflow Pipeline
    ===================================
    Description: The purpose of this pipeline is to automate building a gene database from (DRAM) nucleotide fasta files using DRAM2.

    Usage:
    
    Bring up help menu:
        nextflow run DRAM2.nf --help
    Call genes using input fastas
        nextflow run DRAM2.nf --call --input_fasta_dir <input_fasta_dir> --outdir <outdir> --threads <threads>


    Main DRAM2 Operations:
    --call      : Call genes using prodigal 
    --annotate  : Annotate called genes using user-provided databases
    --distill   : Distill the annotations into an output distillate.xlsx and product.tsv and product.html

    Options:
    --sequence_type     : Type of sequences (nuc or prot).
    --data_type         : Data type of input fasta.
    --input_fasta_dir   : Directory containing input fasta files.
    --outdir            : Output directory path.
    --threads           : Number of threads to use for processing.

    Example:
    nextflow run DRAM2.nf --call --sequence_type nuc --data_type genomic --input_fasta_dir data/input_fasta/ --outdir output --threads 4 
    """.stripIndent()
}

/* Call Help Menu */
def callHelpMessage() {
    log.info """
    DRAM2 Nextflow Pipeline
    ===================================
    Call description: The purpose of DRAM2 --call is to call genes on input FASTA files.

    Usage:

    Call genes using input fastas:
        nextflow run DRAM2.nf --call --input_fasta_dir <input_fasta_dir> --outdir <outdir> --threads <threads>

    Call options:
    --prodigal_mode         'single' or 'meta'
                                Default: 'single'
    --prodigal_tras_table   (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
                            Specify a translation table to use (default: '1').

    Main options:
    --input_fasta           'path/to/fasta/directory/''
                                Directory containing input fasta files.      
                                Default: 'input_fasta/''                  
    --outdir                'path/to/output/directory/''
                                Default: 'DRAM2_output/'
                                Output directory path.
    --threads               Number of threads to use for processing.

    """.stripIndent()
}

/* Annotate Help Menu */
def annotateHelpMessage() {
    log.info """
    DRAM2 Nextflow Pipeline
    ===================================
    Annotate description: The purpose of DRAM2 --annotate is to annotate called genes on input FASTA (faa) files.

    Usage:

    Annotate called genes using input fastas (faa):
        nextflow run DRAM2.nf --annotate --input_fasta_dir <input_fasta_dir> --outdir <outdir> --threads <threads>
    
    Call and Annotate genes using input fastas and KOFAM database:
        nextflow run DRAM2.nf --call --annotate --input_fasta_dir <input_fasta_dir> --outdir <outdir> --threads <threads> --use_kofam

    Annotate options:
    --use_dbset [metabolism_kegg_set|metabolism_set|adjectives|adjectives_kegg]


    --use_[db-name]             [camper, cant_hyd, dbcan, fegenie, kegg, kofam, merops, methyl, heme, pfam, sulfur, uniref]
                                    Specify databases to use. Can use more than one. Can be used in combination with --use_dbset.
    
    --use_dbset                 [metabolism_kegg_set, metabolism_set, adjectives_kegg_set, adjectivs set]
                                    metabolism_kegg_set = kegg, dbcan, merops, pfam, heme
                                    metabolism_set      = kofam, dbcan, merops, pfam, heme
                                    adjectives_kegg_set = kegg, dbcan, merops, pfam, heme, sulfur, camper, methyl, fegenie
                                    adjectives_set      = kofam, dbcan, merops, pfam, heme, sulfur, camper, methyl, fegenie
                                    Only one set can be used. Can be used in combination with --use_[db-name]
    
    --bit_score_threshold       INTEGER
                                    The minimum bit score is calculated by a HMMER or MMseqs search to retain hits.
    --rbh_bit_score_threshold   INTEGER
                                    Minimum bit score of reverse best hits to retain hits.
    --custom_fasta_db PATH      'path/to/fasta/db'
                                    Location of fastas to annotate against, can be used multiple times.
    --custom_hmm_db PATH        'path/to/hmm/db'
                                    Location of HMMs to annotate against, can be used multiple times.
    --custom_hmm_db_cutoffs     'path/to/hmm/cutoffs/db'
                                    Location of file with custom HMM cutoffs and descriptions, can be used multiple times.

    Main options:
    --input_fasta           'path/to/fasta/directory/''
                                Directory containing input fasta files.      
                                Default: 'input_fasta/''                  
    --outdir                'path/to/output/directory/''
                                Default: 'DRAM2_output/'
                                Output directory path.
    --threads               Number of threads to use for processing.

    """.stripIndent()
}

/* Distill Help Menu */
def distillHelpMessage() {
    log.info """
    DRAM2 Nextflow Pipeline
    ===================================
    Distill description: The purpose of DRAM2 --distill is to distill down annotations based on a curated distillation summary form. User's may also provide additional --custom_distillate (TSV forms).

    Usage:

    Distill annotations using input TSV file:
        nextflow run DRAM2.nf --distill --annotations <path/to/annotations.tsv> --outdir <outdir> --threads <threads>
    
    Call and Annotate and Distill genes using input fastas and KOFAM database:
        nextflow run DRAM2.nf --call --annotate --distill --input_fasta_dir <input_fasta_dir> --outdir <outdir> --threads <threads> --use_kofam

    Distill options:
    --annotations_tsv_path PATH     'path/to/annotations.tsv'
                                    This needs to be provided if you are not running distill with --call and --annotate.

    --rrna_path PATH                rRNA output from a dram2 RNA script.
                                        <Description>
    --trna_path PATH                tRNA output from a dram2 annotation.
                                        <Description>
    --show_gene_names               If present, give names of genes instead of counts in genome metabolism summary.
    --custom_distillate             'path/to/custon/distillate.tsv'?
                                        Custom distillate form to add your ownmodules to the metabolism summary. You willneed to read the docs to find the format that this tsv file must take.

    Main options:
    --input_fasta           'path/to/fasta/directory/''
                                Directory containing input fasta files.      
                                Default: 'input_fasta/''                  
    --outdir                'path/to/output/directory/''
                                Default: 'DRAM2_output/'
                                Output directory path.
    --threads               Number of threads to use for processing.

    """.stripIndent()
}

/* Adjectives Help Menu */
def adjectivesHelpMessage() {
    log.info """
    DRAM2 Nextflow Pipeline
    ===================================
    Annotate description: The purpose of DRAM2 --adjectives is to evaluate genes and describe their features.

    Usage:

    Annotate called genes using input fastas (faa):
        nextflow run DRAM2.nf --annotate --input_fasta_dir <input_fasta_dir> --outdir <outdir> --threads <threads>
    
    THESE NEED TO BE RE_DONE
    Adjectives options:
    --annotations_tsv_path PATH     Location of an annotations.tsv. You don't
                                    need to use this option if you are using the
                                    output_dir for dram with a project_config.
                                    If you use this option, you must also use
                                    the force flag to bypass the safeguards that
                                    prevent you from running distill with
                                    insufficient data.
    --adjectives_tsv_path PATH      Location of the output adjectives.tsv. if
                                    you leave this blank the adjectives.tsv file
                                    will be put in the output directory.
    -a, --adjectives TEXT           A list of adjectives, by name, to evaluate.
                                    This limits the number of adjectives that
                                    are evaluated, and is faster.
    -p, --plot_adjectives TEXT      A list of adjectives, by name, to plot. This
                                    limits the number of adjectives that are
                                    plotted and is probably needed for speed.
    -g, --plot_genomes TEXT
    --plot_path PATH                will become a folder of output plots, no
                                    path no plots.
    --strainer_tsv PATH             The path for a tsv that will pass to
                                    strainer to filter genes. The only option at
                                    this time is ‘pgtb’ for positive genes that
                                    are on true bugs.
    --strainer_type PATH            The type of process that should make the
                                    strainer file.
    --debug_ids_by_fasta_to_tsv PATH
                                    This is a tool to debug the list of IDs
                                    found by DRAM it is mostly for experts.
    --user_rules_tsv PATH           This is an optional path to a rules file
                                    with strict formatting. It will overwrite
                                    the original rules file that is stored with
                                    the script.
    --show_rules_path               Show the path to the default rules path.
    --list_name                     List the names for all adjectives_tsv that
                                    are available, you can pass these names to
                                    limit the adjectives that are evaluated
    --list_id                       List the names for all adjectives_tsv that
                                    are available, you can pass these names to
                                    limit the adjectives that are evaluated,
    Main options:
    --input_fasta           'path/to/fasta/directory/''
                                Directory containing input fasta files.      
                                Default: 'input_fasta/''                  
    --outdir                'path/to/output/directory/''
                                Default: 'DRAM2_output/'
                                Output directory path.
    --threads               Number of threads to use for processing.

    """.stripIndent()
}