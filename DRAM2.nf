/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DRAM2: DRAM2 (Distilled and Refined Annotation of Metabolism Version 2) is a tool for annotating metagenomic assembled genomes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Homepage of project
    homePage = https://github.com/WrightonLabCSU/DRAM2

    Author of DRAM2 Nextflow pipeline
    author = Reed Woyda, Rory Flynn
    institution = Colorado State University - Wrighton Lab - Microbial Ecosystems Lab

    Description of project
    description = DRAM2 (Distilled and Refined Annotation of Metabolism Version 2) is a tool for annotating metagenomic assembled genomes.

    Main pipeline script
    mainScript = DRAM2.nf

    version
    v2.0.0
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
include { helpMessage; callHelpMessage; annotateHelpMessage; distillHelpMessage; productHelpMessage; adjectivesHelpMessage; formatKeggHelpMessage; version } from './modules/utils/help.nf'
include { FORMAT_KEGG_DB                                  } from './modules/database/format_kegg_db.nf'

include { RENAME_FASTA                                  } from './modules/call/rename_fasta.nf'
include { CALL_GENES                                    } from './modules/call/call_genes_prodigal.nf'
include { QUAST                                         } from './modules/call/quast.nf'
include { QUAST_COLLECT                                 } from './modules/call/quast_collect.nf'

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
include { ADD_SQL_DESCRIPTIONS as SQL_PFAM              } from './modules/annotate/add_sql_descriptions.nf'

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
include { FEGENIE_HMM_FORMATTER                         } from './modules/annotate/fegenie_hmm_formatter.nf'

include { GENE_LOCS                                     } from './modules/annotate/gene_locs.nf'

include { COMBINE_ANNOTATIONS                           } from './modules/annotate/combine_annotations.nf'
include { COUNT_ANNOTATIONS                             } from './modules/annotate/count_annotations.nf'
include { ADD_ANNOTATIONS                               } from './modules/annotate/add_annotations.nf'
include { MERGE_ANNOTATIONS                             } from './modules/annotate/merge_annotations.nf'
include { GENERATE_GFF_GENBANK                          } from './modules/annotate/generate_gff_genbank.nf'

include { COMBINE_DISTILL                               } from './modules/distill/combine_distill.nf'
include { DISTILL                                       } from './modules/distill/distill.nf'

// This is a placeholder Product process
include { PRODUCT_HEATMAP                               } from './modules/product/product_heatmap.nf'

include { TREES                                         } from './modules/trees/trees.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Help menu and Version menu check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/* Call Help Menu */
if (((params.help || params.h ) && params.call) || params.help_call ){
    callHelpMessage()
    exit 0
}
/* Annotate Help Menu */
else if (((params.help || params.h) && params.call) || params.help_annotate ){
    annotateHelpMessage()
    exit 0
}

/* Distill Help Menu */
else if (((params.help || params.h) && params.distill) || params.help_distill ){
    distillHelpMessage()
    exit 0
}

/* Adjectives Help Menu */
else if (((params.help || params.h) && params.adjectives) || params.help_adjectives ){
    adjectivesHelpMessage()
    exit 0
}

/* Product Help Menu */
else if (((params.help || params.h) && params.product) || params.help_product ){
    productHelpMessage()
    exit 0
}

/* Format KEGG DB Help Menu */
else if (((params.help || params.h) && params.format_kegg) || params.help_format_kegg ){
    formatKeggHelpMessage()
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
def validOptions = ["--call", "--annotate", "--distill_topic", "--distill_ecosystem", "--distill_custom", "--merge_annotations", "--merge_distill", "--rename", "--product", "--format_kegg"]

if ( !params.profile && !params.rename && params.call == 0 && params.annotate == 0 && params.annotations == "" && params.merge_annotations == "" && params.merge_distill == "" && (params.distill_topic == "" || params.distill_ecosystem == "" || params.distill_custom == "" ) && params.format_kegg == 0 ) {
    error("Please provide one of the following options: ${validOptions.join(', ')}")

}

if( !params.profile && params.use_dbset){
    if (!['metabolism_kegg_set', 'metabolism_set', 'adjectives_kegg_set', 'adjectives_set'].contains(params.use_dbset)) {
        error("Invalid parameter '--use_dbset ${params.use_dbset}'. Valid values are 'metabolism_kegg_set', 'metabolism_set', 'adjectives_kegg_set', 'adjectives_set'.")
    }
}

// If --format_kegg and any other options are present, throw an error
if ( params.format_kegg && !(params.call == 0 && params.annotate == 0 && params.annotations == "" && params.merge_annotations == "" && params.merge_distill == "" && (params.distill_topic == "" || params.distill_ecosystem == "" || params.distill_custom == "" ))) {
    error("If you want to format the KEGG database, you must not provide any other options. Format KEGG DB is a standalone process.")
}

if( !params.profile && !params.rename && params.annotations == "" && params.annotate == 0 && (params.distill_topic != "" || params.distill_ecosystem != "" || params.distill_custom != "" )){
    error("If you want to distill, you must provide annotations via --annotations <path/to/file>.")
}

if( params.product && !params.call && !params.annotate && (params.distill_topic == "" || params.distill_ecosystem == "" || params.distill_custom == "" ))
{
    if( params.annotations == "" && params.distillate == "" ){
        error("If you want to generate a product, you must either (1) provide annotations via --annotations <path/to/file> and a distillate --distillate <path/to/file> OR (2) use Call, Annotate and Distill to generate these input files.")
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check for merge_annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Verify the directory exists
if (params.merge_annotations != "") {
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

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Parse DRAM2 ANNOTATE input databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
default_channel = Channel.fromPath(params.distill_dummy_sheet)
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
    ch_fegenie_formatter = file(params.fegenie_hmm_formatter_script)

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

    ch_mmseqs_rbh_script = file(params.mmseqs_rbh_filter_script)

    ch_called_genes_loc_script_faa = file(params.called_genes_loc_script_faa)

    ch_generate_gff_gbk = file(params.ch_generate_gff_gbk_script)

    index_mmseqs = "0"

    if (annotate_kegg == 1) {
        if ( !params.format_kegg && !file(params.kegg_db).exists() ) {
            error("Error: If using --annotate, you must supply prebuilt databases. KEGG database file not found at ${params.kegg_db}")
        }
        else {
            ch_kegg_db = file(params.kegg_db)
        }
        index_mmseqs = "1"
        annotate_list += "${params.kegg_name} "
    }

    if (annotate_kofam == 1) {
        ch_kofam_db = file(params.kofam_db).exists() ? file(params.kofam_db) : error("Error: If using --annotate, you must supply prebuilt databases. KOFAM database file not found at ${params.kofam_db}")
        annotate_list += "${params.kofam_name} "
    }

    if (annotate_dbcan == 1) {
        ch_dbcan_db = file(params.dbcan_db).exists() ? file(params.dbcan_db) : error("Error: If using --annotate, you must supply prebuilt databases. DBCAN database file not found at ${params.dbcan_db}")
        annotate_list += "${params.dbcan_name} "
    }

    if (annotate_camper == 1) {
        ch_camper_hmm_db = file(params.camper_hmm_db).exists() ? file(params.camper_hmm_db) : error("Error: If using --annotate, you must supply prebuilt databases. CAMPER HMM database file not found at ${params.camper_hmm_db}")
        ch_camper_mmseqs_db = file(params.camper_mmseqs_db).exists() ? file(params.camper_mmseqs_db) : error("Error: If using --annotate, you must supply prebuilt databases. CAMPER MMseqs2 database file not found at ${params.camper_mmseqs_db}")
        index_mmseqs = "1"
        annotate_list += "${params.camper_name} "
        ch_camper_mmseqs_list = file(params.camper_mmseqs_list)
    }

    if (annotate_merops == 1) {
        ch_merops_db = file(params.merops_db).exists() ? file(params.merops_db) : error("Error: If using --annotate, you must supply prebuilt databases. MEROPS database file not found at ${params.merops_db}")
        index_mmseqs = "1"
        annotate_list += "${params.merops_name} "
    }

    if (annotate_pfam == 1) {
        ch_pfam_mmseqs_db = file(params.pfam_mmseq_db).exists() ? file(params.pfam_mmseq_db) : error("Error: If using --annotate, you must supply prebuilt databases. PFAM database file not found at ${params.pfam_mmseq_db}")
        index_mmseqs = "1"
        annotate_list += "${params.pfam_name} "
    }

    if (annotate_heme == 1) {
        ch_heme_db = file(params.heme_db).exists() ? file(params.heme_db) : error("Error: If using --annotate, you must supply prebuilt databases. HEME database file not found at ${params.heme_db}")
        annotate_list += "${params.heme_name} "
    }

    if (annotate_sulfur == 1) {
        ch_sulfur_db = file(params.sulfur_db).exists() ? file(params.sulfur_db) : error("Error: If using --annotate, you must supply prebuilt databases. SULURR database file not found at ${params.sulfur_db}")
        annotate_list += "${params.sulfur_name} "
    }

    if (annotate_uniref == 1) {
        ch_uniref_db = file(params.uniref_db).exists() ? file(params.uniref_db) : error("Error: If using --annotate, you must supply prebuilt databases. UNIREF database file not found at ${params.uniref_db}")
        index_mmseqs = "1"
        annotate_list += "${params.uniref_name} "
    }

    if (annotate_methyl == 1) {
        ch_methyl_db = file(params.methyl_db).exists() ? file(params.methyl_db) : error("Error: If using --annotate, you must supply prebuilt databases. METHYL database file not found at ${params.methyl_db}")
        index_mmseqs = "1"
        annotate_list += "${params.methyl_name} "
    }

    if (annotate_fegenie == 1) {
        ch_fegenie_db = file(params.fegenie_db).exists() ? file(params.fegenie_db) : error("Error: If using --annotate, you must supply prebuilt databases. FEGENIE database file not found at ${params.fegenie_db}")
        annotate_list += "${params.fegenie_name} "
    }

    if (annotate_canthyd == 1) {
        ch_canthyd_hmm_db = file(params.canthyd_hmm_db).exists() ? file(params.canthyd_hmm_db) : error("Error: If using --annotate, you must supply prebuilt databases. CANT_HYD HMM database file not found at ${params.canthyd_hmm_db}")
        ch_canthyd_mmseqs_db = file(params.canthyd_mmseqs_db).exists() ? file(params.canthyd_mmseqs_db) : error("Error: If using --annotate, you must supply prebuilt databases. CANT_HYD MMseqs database file not found at ${params.canthyd_mmseqs_db}")
        index_mmseqs = "1"
        annotate_list += "${params.canthyd_name} "
        ch_canthyd_mmseqs_list = file(params.canthyd_mmseqs_list)
    }

    if (annotate_vogdb == 1) {
        ch_vogdb_db = file(params.vog_db).exists() ? file(params.vog_db) : error("Error: If using --annotate, you must supply prebuilt databases. VOG database file not found at ${params.vog_db}")
        annotate_list += "${params.vogdb_name} "
    }

    if (annotate_viral == 1) {
        ch_viral_db = file(params.viral_db).exists() ? file(params.viral_db) : error("Error: If using --annotate, you must supply prebuilt databases. viral database file not found at ${params.viral_db}")
        index_mmseqs = "1"
        annotate_list += "${params.viral_name} "
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

if( params.rename || params.call ){

    ch_generate_gene_locs_script = file(params.generate_gene_locs_script)

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

    if( params.call ){
        // Validate prodigal options for Call
        if (!['single', 'meta'].contains(params.prodigal_mode)) {
            error("Invalid parameter '--prodigal_mode ${params.prodigal_mode}'. Valid values are 'single' or 'meta'.")
        }
        if (!['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25'].contains(params.prodigal_trans_table)) {
            error("Invalid parameter '--prodigal_trans_table ${params.binning_map_mode}'. Valid values are '1','2',...,'25'.")
        }
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
        ch_called_proteins = Channel
            .fromPath(params.input_genes + params.genes_fmt, checkIfExists: true)
            .ifEmpty { exit 1, "If you specify --annotate without --call, you must provide a fasta file of called genes using --input_genes. Cannot find any called gene fasta files matching: ${params.input_genes}\nNB: Path needs to follow pattern: path/to/directory/" }
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
    ch_distill_xlsx_script = file(params.distill_xlsx_script)
    ch_distill_sql_script = file(params.distill_sql_script)

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
            ch_bin_quality = default_channel
        }

        /* Check for input Taxa file */
        if (params.taxa != "") {
            ch_taxa = file(params.taxa).exists() ? file(params.taxa) : error("Error: If using --taxa, you must supply a formatted input file. Taxonomy file not found at ${params.taxa}")
        } else {
            ch_taxa = default_channel
        }

        // Ensure an add_annotations channel is generated if the user specifies --add_annotations
        if( params.add_annotations != ""){
            ch_add_annots = file(params.add_annotations).exists() ? file(params.add_annotations) : error("Error: If using --add_annotations, you must supply a DRAM-formatted annotations file. Annotations file not found at ${params.add_annotations}")
        }

    }

}

if( params.trees ) {

    ch_count_annots_script = file(params.count_annots_script)
    ch_distill_xlsx_script = file(params.distill_xlsx_script)
    ch_distill_sql_script = file(params.distill_sql_script)
    
    //Add in option for --add_trees <list of paths to trees refpkg directories>
    if( params.add_trees ){
        ch_add_trees = file(params.add_trees).exists() ? file(params.add_trees) : error("Error: If using --add_trees, you must supply a path to a directory containing each tree subdirectory. Additional trees directory not found at ${params.add_trees}")
    }
    else{
        ch_add_trees = default_channel
    }    
    


    if( !params.call ){
        if ( params.annotations == "" && params.input_genes == "" ){
            error "If you want to run TREES, you must either use --call to call genes or, provide annotations via --annotations and directory of called genes via --input_genes."
        }
        ch_combined_annotations = Channel
            .fromPath(params.annotations, checkIfExists: true)
            .ifEmpty { exit 1, "If you specify --distill_<topic|ecosystem|custom> without --annotate, you must provide an annotations TSV file (--annotations <path>) with approprite formatting. Cannot find any called gene fasta files matching: ${params.annotations}\nNB: Path needs to follow pattern: path/to/directory/" }

        ch_collected_faa = Channel
            .fromPath(params.input_genes + params.genes_fmt, checkIfExists: true)
            .ifEmpty { exit 1, "If you specify --annotations without --input_genes, with the desire to run trees, you must provide a fasta file of called genes using --input_genes. Cannot find any called gene fasta files matching: ${params.input_genes}\nNB: Path needs to follow pattern: path/to/directory/" }
            .collect()
    }

    ch_tree_data_files = Channel.fromPath(params.tree_data_files)
    ch_trees_scripts = file(params.trees_scripts)

}

if( params.adjectives ){

    ch_adjectives_script = params.adjectives_script

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Create channels for optional topic and/or ecosystem and/or custom sheets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/* Create the default distill topic and ecosystem channels */
def distill_topic_list = ""
def distill_ecosystem_list = ""
def distill_custom_list = ""

distill_default = "0"
distill_carbon = "0"
distill_energy = "0"
distill_misc = "0"
distill_nitrogen = "0"
distill_transport = "0"
distill_camper = "0"

if (params.distill_topic != "" || params.distill_ecosystem != "" || params.distill_custom != "") {
    if (params.distill_topic != "") {
        def validTopics = ['default', 'carbon', 'energy', 'misc', 'nitrogen', 'transport', 'camper']
        def topics = params.distill_topic.split()

        topics.each { topic ->
            if (!validTopics.contains(topic)) {
                error("Invalid distill topic: $topic. Valid values are ${validTopics.join(', ')}")
            }

            // Handle the 'default' case by setting all default topics to "1"
            if (topic == "default") {
                distill_default = "1"
                distill_carbon = "1"
                distill_energy = "1"
                distill_misc = "1"
                distill_nitrogen = "1"
                distill_transport = "1"
                // Ensure that default topics are listed only once
                if (!distill_topic_list.contains("default")) {
                    distill_topic_list += "default (carbon, energy, misc, nitrogen, transport) "
                }
            } else {
                // Handle other topics
                switch (topic) {
                    case "carbon":
                        distill_carbon = "1"
                        if (!distill_topic_list.contains("carbon")) distill_topic_list += "carbon "
                        break
                    case "energy":
                        distill_energy = "1"
                        if (!distill_topic_list.contains("energy")) distill_topic_list += "energy "
                        break
                    case "misc":
                        distill_misc = "1"
                        if (!distill_topic_list.contains("misc")) distill_topic_list += "misc "
                        break
                    case "nitrogen":
                        distill_nitrogen = "1"
                        if (!distill_topic_list.contains("nitrogen")) distill_topic_list += "nitrogen "
                        break
                    case "transport":
                        distill_transport = "1"
                        if (!distill_topic_list.contains("transport")) distill_topic_list += "transport "
                        break
                    case "camper":
                        distill_camper = "1"
                        if (!distill_topic_list.contains("camper")) distill_topic_list += "camper "
                        break
                }
            }
        }

        // Remove the default placeholder if other topics were added in addition to default
        if (distill_default == "1" && topics.size() > 1) {
            distill_topic_list = distill_topic_list.replace("default (carbon, energy, misc, nitrogen, transport) ", "")
            if (!distill_topic_list.contains("carbon")) distill_topic_list = "carbon, energy, misc, nitrogen, transport, " + distill_topic_list.trim()
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

    if (distill_camper == "1") {
        ch_distill_camper = file(params.distill_camper_sheet).exists() ? file(params.distill_camper_sheet) : error("Error: If using --distill_topic camper, you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

    } else{
        ch_distill_camper = default_channel
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
    /*
    if (params.distill_custom != "") {
        ch_distill_custom = file(params.distill_custom).exists() ? file(params.distill_custom) : error("Error: If using --distill_custom <path/to/TSV>, you must have the preformatted custom distill sheet in the provided file: ${params.distill_custom}.")
        distill_custom_list = "params.distill_custom"
    }else{
        ch_distill_custom = default_channel
    }
    */
    if (params.distill_custom != "") {
    // Verify the directory exists
    def custom_distill_dir = file(params.distill_custom)
    if (!custom_distill_dir.exists()) {
        error "Error: The specified directory for merging annotations (--distill_custom) does not exist: ${params.distill_custom}"
    }
    else{
        ch_distill_custom_collected = default_channel
    }

    // Verify the directory contains .tsv files
    def tsv_files = custom_distill_dir.list().findAll { it.endsWith('.tsv') }
    if (tsv_files.isEmpty()) {
        error "Error: The specified directory for merging annotations (--distill_custom) does not contain any .tsv files: ${params.distill_custom}"
    }
    else{
            // Create a channel with the paths to the .tsv files
    Channel
        .from(tsv_files.collect { custom_distill_dir.toString() + '/' + it })
        .set { ch_distill_custom }
    Channel.empty()
        .mix( ch_distill_custom )
        .collect()
        .set { ch_distill_custom_collected }

    }
    }
    else{
        ch_distill_custom_collected = default_channel
    }

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Create channels for Formatting Kegg DB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
if ( params.format_kegg ) {
    ch_format_kegg_db_script = file(params.format_kegg_db_script)
    ch_kegg_pep = file(params.kegg_pep_loc).exists() ? file(params.kegg_pep_loc) : error("Error: when using running format_kegg, the kegg_pep_loc file must exists. KEGG pep file not found at ${params.kegg_pep_loc}")
    /*
    if ( params.gene_ko_link_loc != "" ) {
        ch_gene_ko_link = file(params.gene_ko_link_loc).exists() ? file(params.gene_ko_link_loc) : error("Error: If supplying gene_ko_link_loc, the file must exists. Gene-KO link file not found at ${params.gene_ko_link_loc}")
    }
    else {
        ch_gene_ko_link = ['']
    
    }
    */
    
    // ch_gene_ko_link = Channel.fromPath(params.gene_ko_link_loc).ifEmpty('NAN FILE')
    
    if ( params.gene_ko_link_loc != "" ) {
        ch_gene_ko_link = Channel.fromPath(params.gene_ko_link_loc).ifEmpty('')
    }
    else {
       ch_gene_ko_link = Channel.empty().ifEmpty('')
    }
    if ( params.annotate || annotate_kegg != 1 ) {
        ch_kegg_db = file(params.kegg_db)
    }
    if ( params.kegg_download_date ) {
        kegg_download_date = params.kegg_download_date
    }
    else {
        kegg_download_date = "''"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Print run info to command line
        Various run info for various DRAM2 options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// This is just a catch-all for now - NEED to generate others for various options
if( params.annotate && params.call == "" && (params.distill_ecosystem =="" || params.distill_custom =="" || params.distill_topic =="" )){
    log.info """
            DRAM2 Nextflow
            ===================================
            fastas       : ${params.input_fasta}
            outdir       : ${params.outdir}
            threads      : ${params.threads}
            annotate     : ${params.annotate ? 'true' : 'false'}
            databases    : ${annotate_list}

            """
            .stripIndent()
}else if( params.call && params.annotate && (params.distill_ecosystem !="" || params.distill_custom !="" || params.distill_topic !="" )){
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
            databases    : ${annotate_list}
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


    /* If we are formatting kegg, we do that and then exit the program */
    if ( params.format_kegg ) {
        FORMAT_KEGG_DB( ch_kegg_pep, ch_gene_ko_link, ch_format_kegg_db_script, kegg_download_date )
        return
    }

    /* Rename fasta headers
        Process 1-by-1
    */

    if( params.rename ) {
        RENAME_FASTA( ch_input_fastas )
        ch_fasta = RENAME_FASTA.out.renamed_fasta
    }
    else if (params.call) {
        ch_fasta = ch_input_fastas
    }

    if( params.call ){
        // Call genes using Prodigal on the input fasta file(s) 1-by-1
        CALL_GENES ( ch_fasta, ch_generate_gene_locs_script )
        ch_called_genes = CALL_GENES.out.prodigal_fna
        ch_called_proteins = CALL_GENES.out.prodigal_faa
        ch_gene_locs = CALL_GENES.out.prodigal_locs_tsv
        ch_gene_gff = CALL_GENES.out.prodigal_gff
        ch_filtered_fasta = CALL_GENES.out.prodigal_filtered_fasta

        // Collect all individual fasta to pass to quast
        ch_called_proteins
            .map { tuple -> tuple[1] }  // Extract only the file path from each tuple
            .collect()                  // Collect all paths into a list
            .set { ch_collected_faa }   // Set the resulting list to ch_collected_faa

        // Collect all individual fasta to pass to quast
        Channel.empty()
            .mix( ch_called_genes  )
            .collect()
            .set { ch_collected_fna }

        // Collect all individual fasta to pass to quast
        Channel.empty()
            .mix( ch_filtered_fasta, ch_gene_gff  )
            .collect()
            .set { ch_collected_fasta }

        // Run QUAST on individual FASTA file combined with respective GFF
        QUAST( ch_collected_fasta )
        ch_quast_stats = QUAST.out.quast_collected_out

        // Run tRNAscan-SE on each fasta to identify tRNAs
        TRNA_SCAN( ch_fasta )
        ch_trna_scan = TRNA_SCAN.out.trna_scan_out
        // Collect all sample formatted tRNA files
        Channel.empty()
            .mix( ch_trna_scan )
            .collect()
            .set { ch_collected_tRNAs }

        // Run TRNA_COLLECT to generate a combined TSV for all fastas
        TRNA_COLLECT( ch_collected_tRNAs )
        ch_trna_sheet = TRNA_COLLECT.out.trna_collected_out

        // Run barrnap on each fasta to identify rRNAs
        RRNA_SCAN( ch_fasta )
        ch_rrna_scan = RRNA_SCAN.out.rrna_scan_out
        Channel.empty()
            .mix( ch_rrna_scan )
            .collect()
            .set { ch_collected_rRNAs }

        // Run RRNA_COLLECT to generate a combined TSV for all fastas
        RRNA_COLLECT( ch_collected_rRNAs )
        ch_rrna_sheet = RRNA_COLLECT.out.rrna_collected_out
        ch_rrna_combined = RRNA_COLLECT.out.rrna_combined_out


    }
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Merge Annotations
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if( params.merge_annotations != "" ){
        MERGE_ANNOTATIONS( ch_merge_annotations_collected )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Annotation
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if( params.annotate ){

        // Define empty channel to populate with annotation results
        def formattedOutputChannels = channel.of()

        if( params.call == 0){
            GENE_LOCS( ch_called_proteins, ch_called_genes_loc_script_faa )
            ch_gene_locs = GENE_LOCS.out.prodigal_locs_tsv
        }

        // Here we will create mmseqs2 index files for each of the inputs if we are going to do a mmseqs2 database
        if( index_mmseqs == "1" ){
            // Use MMSEQS2 to index each called genes protein file
            MMSEQS_INDEX( ch_called_proteins )
            ch_mmseqs_query = MMSEQS_INDEX.out.mmseqs_index_out
        }
        // KEGG annotation
        if( annotate_kegg == 1 ){
            ch_combined_query_locs_kegg = ch_mmseqs_query.join(ch_gene_locs)
            MMSEQS_SEARCH_KEGG( ch_combined_query_locs_kegg, ch_kegg_db, params.bit_score_threshold, params.rbh_bit_score_threshold, ch_dummy_sheet, params.kegg_name, ch_mmseqs_script, ch_mmseqs_rbh_script )
            ch_kegg_unformatted = MMSEQS_SEARCH_KEGG.out.mmseqs_search_formatted_out

            SQL_KEGG(ch_kegg_unformatted, params.kegg_name, ch_sql_descriptions_db, ch_sql_parser)
            ch_kegg_formatted = SQL_KEGG.out.sql_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_kegg_formatted)
        }
        // KOFAM annotation
        if( annotate_kofam == 1 ){
            HMM_SEARCH_KOFAM ( ch_called_proteins, params.kofam_e_value, ch_kofam_db )
            ch_kofam_hmms = HMM_SEARCH_KOFAM.out.hmm_search_out

            PARSE_HMM_KOFAM ( ch_kofam_hmms, ch_parse_hmmsearch )
            ch_kofam_parsed = PARSE_HMM_KOFAM.out.parsed_hmm

            ch_combined_hits_locs_kofam = ch_kofam_parsed.join(ch_gene_locs)
            KOFAM_HMM_FORMATTER ( ch_combined_hits_locs_kofam, params.kofam_top_hit, ch_kofam_list, ch_kofam_formatter )
            ch_kofam_formatted = KOFAM_HMM_FORMATTER.out.kofam_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_kofam_formatted)
        }
        // PFAM annotation
        if( annotate_pfam == 1 ){
            ch_combined_query_locs_pfam = ch_mmseqs_query.join(ch_gene_locs)
            MMSEQS_SEARCH_PFAM( ch_combined_query_locs_pfam, ch_pfam_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, ch_dummy_sheet, params.pfam_name, ch_mmseqs_script, ch_mmseqs_rbh_script )
            ch_pfam_unformatted = MMSEQS_SEARCH_PFAM.out.mmseqs_search_formatted_out

            SQL_PFAM(ch_pfam_unformatted, params.pfam_name, ch_sql_descriptions_db, ch_sql_parser)
            ch_pfam_formatted = SQL_PFAM.out.sql_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_pfam_formatted)
        }
        // dbCAN annotation
        if( annotate_dbcan == 1 ){

            HMM_SEARCH_DBCAN ( ch_called_proteins, params.dbcan_e_value , ch_dbcan_db)
            ch_dbcan_hmms = HMM_SEARCH_DBCAN.out.hmm_search_out

            PARSE_HMM_DBCAN ( ch_dbcan_hmms, ch_parse_hmmsearch )
            ch_dbcan_parsed = PARSE_HMM_DBCAN.out.parsed_hmm

            ch_combined_hits_locs_dbcan = ch_dbcan_parsed.join(ch_gene_locs)
            DBCAN_HMM_FORMATTER ( ch_combined_hits_locs_dbcan, params.dbcan_top_hit, params.dbcan_name, ch_dbcan_formatter, ch_sql_parser, ch_sql_descriptions_db )
            ch_dbcan_formatted = DBCAN_HMM_FORMATTER.out.sql_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_dbcan_formatted)
        }
        // CAMPER annotation
        if (annotate_camper == 1){
            // HMM
            HMM_SEARCH_CAMPER ( ch_called_proteins, params.camper_e_value , ch_camper_hmm_db)
            ch_camper_hmms = HMM_SEARCH_CAMPER.out.hmm_search_out

            PARSE_HMM_CAMPER ( ch_camper_hmms, ch_parse_hmmsearch )
            ch_camper_parsed = PARSE_HMM_CAMPER.out.parsed_hmm

            ch_combined_hits_locs_camper = ch_camper_parsed.join(ch_gene_locs)
            CAMPER_HMM_FORMATTER ( ch_combined_hits_locs_camper, params.camper_top_hit, ch_camper_hmm_list, ch_camper_formatter )
            ch_camper_hmm_formatted = CAMPER_HMM_FORMATTER.out.camper_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_camper_hmm_formatted)

            // MMseqs
            ch_combined_query_locs_camper = ch_mmseqs_query.join(ch_gene_locs)
            MMSEQS_SEARCH_CAMPER( ch_combined_query_locs_camper, ch_camper_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, ch_camper_mmseqs_list, params.camper_name, ch_mmseqs_script, ch_mmseqs_rbh_script )
            ch_camper_mmseqs_formatted = MMSEQS_SEARCH_CAMPER.out.mmseqs_search_formatted_out

            formattedOutputChannels = formattedOutputChannels.mix(ch_camper_mmseqs_formatted)
        }
        // FeGenie annotation
        if (annotate_fegenie == 1){
            HMM_SEARCH_FEGENIE ( ch_called_proteins,  params.fegenie_e_value, ch_fegenie_db )
            ch_fegenie_hmms = HMM_SEARCH_FEGENIE.out.hmm_search_out

            PARSE_HMM_FEGENIE ( ch_fegenie_hmms, ch_parse_hmmsearch )
            ch_fegenie_parsed = PARSE_HMM_FEGENIE.out.parsed_hmm

            ch_combined_hits_locs_fegenie = ch_fegenie_parsed.join(ch_gene_locs)
            FEGENIE_HMM_FORMATTER ( ch_combined_hits_locs_fegenie, ch_fegenie_formatter )
            ch_fegenie_formatted = FEGENIE_HMM_FORMATTER.out.fegenie_formatted_hits
            formattedOutputChannels = formattedOutputChannels.mix(ch_fegenie_formatted)
        }
        // Methyl annotation
        if (annotate_methyl == 1){
            ch_combined_query_locs_methyl = ch_mmseqs_query.join(ch_gene_locs)
            MMSEQS_SEARCH_METHYL( ch_combined_query_locs_methyl, ch_methyl_db, params.bit_score_threshold, params.rbh_bit_score_threshold, ch_dummy_sheet, params.methyl_name, ch_mmseqs_script, ch_mmseqs_rbh_script )
            ch_methyl_mmseqs_formatted = MMSEQS_SEARCH_METHYL.out.mmseqs_search_formatted_out

            formattedOutputChannels = formattedOutputChannels.mix(ch_methyl_mmseqs_formatted)
        }
        // CANT-HYD annotation
        if (annotate_canthyd == 1){
            // MMseqs
            ch_combined_query_locs_canthyd = ch_mmseqs_query.join(ch_gene_locs)
            MMSEQS_SEARCH_CANTHYD( ch_combined_query_locs_canthyd, ch_canthyd_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, ch_canthyd_mmseqs_list, params.canthyd_name, ch_mmseqs_script, ch_mmseqs_rbh_script )
            ch_canthyd_mmseqs_formatted = MMSEQS_SEARCH_CANTHYD.out.mmseqs_search_formatted_out

            formattedOutputChannels = formattedOutputChannels.mix(ch_canthyd_mmseqs_formatted)

            //HMM
            HMM_SEARCH_CANTHYD ( ch_called_proteins, params.canthyd_e_value , ch_canthyd_hmm_db)
            ch_canthyd_hmms = HMM_SEARCH_CANTHYD.out.hmm_search_out

            PARSE_HMM_CANTHYD ( ch_canthyd_hmms, ch_parse_hmmsearch )
            ch_canthyd_parsed = PARSE_HMM_CANTHYD.out.parsed_hmm

            ch_combined_hits_locs_canthyd = ch_canthyd_parsed.join(ch_gene_locs)
            CANTHYD_HMM_FORMATTER ( ch_combined_hits_locs_canthyd, params.canthyd_top_hit, ch_canthyd_hmm_list, ch_canthyd_formatter )
            ch_canthyd_hmm_formatted = CANTHYD_HMM_FORMATTER.out.canthyd_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_canthyd_hmm_formatted)

        }
        // Sulfur annotation
        if (annotate_sulfur == 1){
            HMM_SEARCH_SULFUR ( ch_called_proteins,  params.sulfur_e_value, ch_sulfur_db )
            ch_sulfur_hmms = HMM_SEARCH_SULFUR.out.hmm_search_out

            PARSE_HMM_SULFUR ( ch_sulfur_hmms, ch_parse_hmmsearch )
            ch_sulfur_parsed = PARSE_HMM_SULFUR.out.parsed_hmm

            ch_combined_hits_locs_sulfur = ch_sulfur_parsed.join(ch_gene_locs)
            SULFUR_HMM_FORMATTER ( ch_combined_hits_locs_sulfur, ch_sulfur_formatter )
            ch_sulfur_formatted = SULFUR_HMM_FORMATTER.out.sulfur_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_sulfur_formatted)
        }
        // MEROPS annotation
        if (annotate_merops == 1){
            ch_combined_query_locs_merops = ch_mmseqs_query.join(ch_gene_locs)
            MMSEQS_SEARCH_MEROPS( ch_combined_query_locs_merops, ch_merops_db, params.bit_score_threshold, params.rbh_bit_score_threshold, ch_dummy_sheet, params.merops_name, ch_mmseqs_script, ch_mmseqs_rbh_script )
            ch_merops_unformatted = MMSEQS_SEARCH_MEROPS.out.mmseqs_search_formatted_out

            SQL_MEROPS(ch_merops_unformatted, params.merops_name, ch_sql_descriptions_db, ch_sql_parser)
            ch_merops_formatted = SQL_MEROPS.out.sql_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_merops_formatted)
        }
        // Uniref annotation
        if (annotate_uniref == 1){
            ch_combined_query_locs_uniref = ch_mmseqs_query.join(ch_gene_locs)
            MMSEQS_SEARCH_UNIREF( ch_combined_query_locs_uniref, ch_uniref_db, params.bit_score_threshold, params.rbh_bit_score_threshold, ch_dummy_sheet, params.uniref_name, ch_mmseqs_script, ch_mmseqs_rbh_script )
            ch_uniref_unformatted = MMSEQS_SEARCH_UNIREF.out.mmseqs_search_formatted_out

            SQL_UNIREF(ch_uniref_unformatted, params.uniref_name, ch_sql_descriptions_db, ch_sql_parser)
            ch_uniref_formatted = SQL_UNIREF.out.sql_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_uniref_formatted)
        }
        // VOGdb annotation
        if (annotate_vogdb == 1){
            HMM_SEARCH_VOG ( ch_called_proteins, params.vog_e_value , ch_vogdb_db )
            ch_vog_hmms = HMM_SEARCH_VOG.out.hmm_search_out

            PARSE_HMM_VOG ( ch_vog_hmms, ch_parse_hmmsearch )
            ch_vog_parsed = PARSE_HMM_VOG.out.parsed_hmm

            ch_combined_hits_locs_vog = ch_vog_parsed.join(ch_gene_locs)
            VOG_HMM_FORMATTER ( ch_combined_hits_locs_vog, params.vog_top_hit, params.vogdb_name, ch_vog_formatter, ch_sql_parser, ch_sql_descriptions_db )
            ch_vog_formatted = VOG_HMM_FORMATTER.out.vog_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_vog_formatted)
        }
        // Viral annotation
        if (annotate_viral == 1){
            ch_combined_query_locs_viral = ch_mmseqs_query.join(ch_gene_locs)
            MMSEQS_SEARCH_VIRAL( ch_combined_query_locs_viral, ch_viral_db, params.bit_score_threshold,  params.rbh_bit_score_threshold,ch_dummy_sheet, params.viral_name, ch_mmseqs_script, ch_mmseqs_rbh_script )
            ch_viral_unformatted = MMSEQS_SEARCH_VIRAL.out.mmseqs_search_formatted_out

            SQL_VIRAL(ch_viral_unformatted, params.viral_name, ch_sql_descriptions_db, ch_sql_parser)
            ch_viral_formatted = SQL_VIRAL.out.sql_formatted_hits

            formattedOutputChannels = formattedOutputChannels.mix(ch_viral_formatted)
        }

        // Collect all formatted annotation output files
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
            ADD_ANNOTATIONS( ch_updated_taxa_annots, ch_add_annots )
            ch_final_annots = ADD_ANNOTATIONS.out.combined_annots_out

            // If the user wants to run trees, do it before we count the annotations
            if( params.trees ){
                TREES( ch_final_annots, params.trees_list, ch_collected_faa, ch_tree_data_files, ch_trees_scripts, ch_add_trees )
            }

            COUNT_ANNOTATIONS ( ch_final_annots, ch_count_annots_script, ch_distill_sql_script  )
            ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts
            ch_annotations_sqlite3 = COUNT_ANNOTATIONS.out.annotations_sqlite3
        }
        else{
            // If the user wants to run trees, do it before we count the annotations
            if( params.trees ){
                TREES( ch_updated_taxa_annots, params.trees_list, ch_collected_faa, ch_tree_data_files, ch_trees_scripts, ch_add_trees )
            }

            ch_final_annots = ch_updated_taxa_annots
            COUNT_ANNOTATIONS ( ch_final_annots, ch_count_annots_script, ch_distill_sql_script  )
            ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts
            ch_annotations_sqlite3 = COUNT_ANNOTATIONS.out.annotations_sqlite3
        }

        if( params.generate_gff || params.generate_gbk ){
            GENERATE_GFF_GENBANK( ch_collected_fna, params.database_list, ch_final_annots, ch_generate_gff_gbk )
        }


    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Solo Trees
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if(params.annotations != "" && params.input_genes != "" && params.distill_topic == "" && params.distill_ecosystem == "" && params.distill_custom == ""){
        
        if( params.trees ){
            TREES( ch_combined_annotations, params.trees_list, ch_collected_faa, ch_tree_data_files, ch_trees_scripts, ch_add_trees )
            ch_combined_annotations = TREES.out.updated_annotations

            COUNT_ANNOTATIONS ( ch_combined_annotations, ch_count_annots_script, ch_distill_sql_script  )
            ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts
            ch_annotations_sqlite3 = COUNT_ANNOTATIONS.out.annotations_sqlite3
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
        if( params.call == 0 ){
            if( params.trnas != "" ){
                TRNA_COLLECT( ch_collected_tRNAs )
                ch_trna_sheet = TRNA_COLLECT.out.trna_collected_out
            }else{
                ch_trna_sheet = default_channel
            }
            if( params.rrnas != "" ){
            RRNA_COLLECT( ch_collected_rRNAs )
                ch_rrna_sheet = RRNA_COLLECT.out.rrna_collected_out
                ch_rrna_combined = RRNA_COLLECT.out.rrna_combined_out
            }else{
                ch_rrna_sheet = default_channel
                ch_rrna_combined = default_channel
            }
            ch_quast_stats = default_channel
        }
        // If the user did not annotate and provided taxa and/or bin quality, add it to annotations.
        if( params.annotate == 0 ){
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
            if( params.add_annotations != "" ){
                // Add additional annotations if user provided them
                ADD_ANNOTATIONS( ch_updated_taxa_annots, ch_add_annots )
                ch_final_annots = ADD_ANNOTATIONS.out.combined_annots_out

                if( params.trees ){
                    TREES( ch_combined_annotations, params.trees_list, ch_collected_faa, ch_tree_data_files, ch_trees_scripts, ch_add_trees )
                    ch_trees_updated_annots = TREES.out.updated_annotations
                }
                else{
                    ch_trees_updated_annots = ch_final_annots
                }
                ch_final_annots = ch_trees_updated_annots
                // Count annotations per sample
                COUNT_ANNOTATIONS ( ch_final_annots, ch_count_annots_script, ch_distill_sql_script )
                ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts
                ch_annotations_sqlite3 = COUNT_ANNOTATIONS.out.annotations_sqlite3
            }
            else{
                if( params.trees ){
                    ch_collected_faa.view()
                    TREES( ch_updated_taxa_annots, params.trees_list, ch_collected_faa, ch_tree_data_files, ch_trees_scripts, ch_add_trees )
                    ch_trees_updated_annots = TREES.out.updated_annotations
                }
                else{
                    ch_trees_updated_annots = ch_updated_taxa_annots
                }
                ch_final_annots = ch_trees_updated_annots
                // Count annotations per sample
                COUNT_ANNOTATIONS ( ch_final_annots, ch_count_annots_script, ch_distill_sql_script )
                ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts
                ch_annotations_sqlite3 = COUNT_ANNOTATIONS.out.annotations_sqlite3
            }ne

        }

        // Combine the individual user-specified distill sheets into a single channel
        COMBINE_DISTILL(ch_distill_carbon, ch_distill_energy, ch_distill_misc, ch_distill_nitrogen, ch_distill_transport, ch_distill_ag, ch_distill_eng_sys, ch_distill_camper, ch_distill_custom_collected )
        ch_combined_distill_sheets = COMBINE_DISTILL.out.ch_combined_distill_sheets

        // Generate multi-sheet XLSX document containing annotations included in user-specified distillate speadsheets

        DISTILL( ch_final_annots, ch_combined_distill_sheets, ch_annotation_counts, ch_quast_stats, ch_rrna_sheet, ch_rrna_combined, ch_trna_sheet, ch_distill_xlsx_script, ch_annotations_sqlite3 )
        ch_distillate = DISTILL.out.distillate
    }
    else if (params.product) {  // if not running any distillates but still running product, just pass the cli annotations file to the next step
        ch_final_annots = file(params.annotations)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Product
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if( params.product ){

        PRODUCT_HEATMAP( ch_final_annots, params.groupby_column )

    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Adjectives
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    /*
    if( params.adjectives ){
        ADJECTIVES( ch_final_annots, ch_adjectives_script )

    }
    */
}
