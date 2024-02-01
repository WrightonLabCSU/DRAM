/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DRAM2: <Description>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Homepage of project 
    homePage = <GitHub homepage>
    
    Author of DRAM2 Nextflow pipeline
    author = Reed Woyda, Rory Flynn
    institutioon = Colorado State University - Wrighton Lab

    Description of project
    description = <Description>

    Main pipeline script
    mainScript = DRAM2.nf
    
    version
    v0.0.1
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
    Help menu and Version menu check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/* Call Help Menu */
if (((params.help) || (params.h)) && params.call ){
    callHelpMessage()
    exit 0
}
/* Annotate Help Menu */
else if (((params.help) || (params.h)) && params.annotate ){
    annotateHelpMessage()
    exit 0
}

/* Distill Help Menu */
else if (((params.help) || (params.h)) && params.distill ){
    distillHelpMessage()
    exit 0
}

/* Adjectives Help Menu */
else if (((params.help) || (params.h)) && params.adjectives ){
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

if (!params.call && !params.annotate && (params.distill_topic == "" || params.distill_ecosystem == "" || params.distill_custom == "" )) {
    error("Please provide one of the following options: ${validOptions.join(', ')}")
}

if( params.use_dbset){
    if (!['metabolism_kegg_set', 'metabolism_set', 'adjectives_kegg_set', 'adjectives_set'].contains(params.use_dbset)) {
        error("Invalid parameter '--use_dbset ${params.use_dbset}'. Valid values are 'metabolism_kegg_set', 'metabolism_set', 'adjectives_kegg_set', 'adjectives_set'.")
    }
}

//Add in other checks

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

include { MMSEQS2                                       } from './modules/mmseqs2.nf'

include { HMM_SEARCH as HMM_SEARCH_KOFAM                } from './modules/annotate/hmmsearch.nf'
include { PARSE_HMM as PARSE_HMM_KOFAM                  } from './modules/annotate/parse_hmmsearch.nf'

include { HMM_SEARCH as HMM_SEARCH_DBCAN                } from './modules/annotate/hmmsearch.nf'
include { PARSE_HMM as PARSE_HMM_DBCAN                  } from './modules/annotate/parse_hmmsearch.nf'

include { GENERIC_HMM_FORMATTER                         } from './modules/annotate/generic_hmm_formatter.nf'
include { KEGG_HMM_FORMATTER                            } from './modules/annotate/kegg_hmm_formatter.nf'
include { KOFAM_HMM_FORMATTER                           } from './modules/annotate/kofam_hmm_formatter.nf'
include { DBCAN_HMM_FORMATTER                           } from './modules/annotate/dbcan_hmm_formatter.nf'

include { INDEX as KEGG_INDEX                           } from './modules/index.nf'

include { COMBINE_ANNOTATIONS                           } from './modules/annotate/combine_annotations.nf'
include { COUNT_ANNOTATIONS                             } from './modules/annotate/count_annotations.nf'

include { COMBINE_DISTILL                               } from './modules/distill/combine_distill.nf'
include { DISTILL_SUMMARY                               } from './modules/distill/distill_summary.nf'
include { DISTILL_FINAL                                 } from './modules/distill/distill_final.nf'

include { PRODUCT_HEATMAP                               } from './modules/product/product_heatmap.nf'


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
annotate_cant_hyd = 0
annotate_heme = 0
annotate_sulfur = 0
annotate_methyl = 0
annotate_merops = 0

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
if( params.use_cant_hyd ){
    annotate_cant_hyd = 1
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

if( params.annotate ){
    //This is just temporary - want these in the containers eventually
    ch_parse_hmmsearch = file(params.parse_hmmsearch_script)
    ch_kegg_formatter = file(params.kegg_formatter_script)
    ch_combine_annot_script = file(params.combine_annotations_script)
    ch_count_annots_script = file(params.count_annots_script)
    ch_distill_summary_script = file(params.distill_summary_script)
    ch_distill_final_script = file(params.distill_final_script)
    ch_kofam_formatter = file(params.kofam_hmm_formatter_script)
    ch_kofam_list = file(params.kofam_list)
    ch_dbcan_formatter = file(params.dbcan_hmm_formatter_script)
    ch_dbcan_fam = file(params.dbcan_fam_activities)
    ch_dbcan_subfam = file(params.dbcan_subfam_activities)

    if (annotate_kegg == 1) {
        ch_kegg_db = file(params.kegg_db).exists() ? file(params.kegg_db) : error("Error: If using --annotate, you must supply prebuilt databases. KEGG database file not found at ${params.kegg_db}")
    } else {
        ch_kegg_db = []
    }

    if (annotate_kofam == 1) {
        ch_kofam_db = file(params.kofam_db).exists() ? file(params.kofam_db) : error("Error: If using --annotate, you must supply prebuilt databases. KOFAM database file not found at ${params.kofam_db}")
    } else {
        ch_kofam_db = []
    }

    if (annotate_dbcan == 1) {
        ch_dbcan_db = file(params.dbcan_db).exists() ? file(params.dbcan_db) : error("Error: If using --annotate, you must supply prebuilt databases. DBCAN database file not found at ${params.dbcan_db}")
    } else {
        ch_dbcan_db = []
    }

    if (annotate_camper == 1) {
        ch_camper_db = file(params.camper_db).exists() ? file(params.camper_db) : error("Error: If using --annotate, you must supply prebuilt databases. CAMPER database file not found at ${params.camper_db}")
    } else {
        ch_camper_db = []
    }

    if (annotate_merops == 1) {
        ch_merops_db = file(params.merops_db).exists() ? file(params.merops_db) : error("Error: If using --annotate, you must supply prebuilt databases. MEROPS database file not found at ${params.merops_db}")
    } else {
        ch_merops_db = []
    }

    if (annotate_pfam == 1) {
        ch_pfam_db = file(params.pfam_db).exists() ? file(params.pfam_db) : error("Error: If using --annotate, you must supply prebuilt databases. PFAM database file not found at ${params.pfam_db}")
    } else {
        ch_pfam_db = []
    }

    if (annotate_heme == 1) {
        ch_heme_db = file(params.heme_db).exists() ? file(params.heme_db) : error("Error: If using --annotate, you must supply prebuilt databases. HEME database file not found at ${params.heme_db}")
    } else {
        ch_heme_db = []
    }

    if (annotate_sulfur == 1) {
        ch_sulfur_db = file(params.sulfur_db).exists() ? file(params.sulfur_db) : error("Error: If using --annotate, you must supply prebuilt databases. SULFUR database file not found at ${params.sulfur_db}")
    } else {
        ch_sulfur_db = []
    }

    if (annotate_methyl == 1) {
        ch_methyl_db = file(params.methyl_db).exists() ? file(params.methyl_db) : error("Error: If using --annotate, you must supply prebuilt databases. METHYL database file not found at ${params.methyl_db}")
    } else {
        ch_methyl_db = []
    }

    if (annotate_fegenie == 1) {
        ch_fegenie_db = file(params.fegenie_db).exists() ? file(params.fegenie_db) : error("Error: If using --annotate, you must supply prebuilt databases. FEGENIE database file not found at ${params.fegenie_db}")
    } else {
        ch_fegenie_db = []
    }

    if (annotate_cant_hyd == 1) {
        ch_cant_hyd_db = file(params.cant_hyd_db).exists() ? file(params.cant_hyd_db) : error("Error: If using --annotate, you must supply prebuilt databases. CANT_HYD database file not found at ${params.cant_hyd_db}")
    } else {
        ch_cant_hyd_db = []
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
// Here we will check for the various kinds of inputs
// The structure of this will change
// For example, we will be able to ingest only metaT data and call genes on that
// For example, users may use a csv or put a path to a directory
if( params.call ){
    // If calling genes, then create a channel called ch_fastas.
    if ( params.input_fasta ) {
        ch_fastas = Channel
            .fromPath(params.input_fasta + params.fasta_fmt, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find any fasta files matching: ${params.input_fasta}\nNB: Path needs to follow pattern: path/to/directory/" }
    } else {
        ch_fastas = Channel
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

    // Convert ch_fastas into a tuple: [samplename, file]
    ch_fastas.map {
        sampleName = it.getName().replaceAll(/\.[^.]+$/, '').replaceAll(/\./, '-')
        tuple(sampleName, it)
    }.set{fastas}
}

if( params.annotate ){

    if( params.call == 0 ){
        ch_called_genes = Channel
            .fromPath(params.input_genes + params.genes_fmt, checkIfExists: true)
            .ifEmpty { exit 1, "If you specify --annotate without --call, you must provide a fasta files of called genes. Cannot find any called gene fasta files matching: ${params.input_genes}\nNB: Path needs to follow pattern: path/to/directory/" }
    }
    // Convert the input_genes into a tuple: [samplename, file]
    if( params.input_genes != 0 ){
        called_proteins = Channel.fromPath(ch_called_genes).map {
            sampleName = it.getName().replaceAll(/\.[^.]+$/, '').replaceAll(/\./, '-')
            tuple(sampleName, it)
        }
        called_proteins = ch_called_genes
    }
}
/*
if( params.merge ){
    
}
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Create channels for optional topic and/or ecosystem and/or custom sheets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/* Create the default distill topic and ecosystem channels */
default_channel = Channel.fromPath(params.distill_dummy_sheet)

if (params.distill_topic != "" || params.distill_ecosystem != "" || params.distill_custom != "") {    
    if (params.distill_topic != "") {
        distill_default = 0
        distill_carbon = 0
        distill_energy = 0
        distill_misc = 0
        distill_nitrogen = 0
        distill_transport = 0

        def validTopics = ['default', 'carbon', 'energy', 'misc', 'nitrogen', 'transport']
        def topics = params.distill_topic.split()

        // Create a list to store the generated channels
        def topicChannels = []

        topics.each { topic ->
            if (!validTopics.contains(topic)) {
                error("Invalid distill topic: $topic. Valid values are ${validTopics.join(', ')}")
            }

            switch (topic) {
                case "default":
                    distill_carbon = 1
                    distill_energy = 1
                    distill_misc = 1
                    distill_nitrogen = 1
                    distill_transport = 1
                    break
                case "carbon":
                    distill_carbon = 1
                    break
                case "energy":
                    distill_energy = 1
                    break
                case "misc":
                    distill_misc = 1
                    break
                case "nitrogen":
                    distill_nitrogen = 1
                    break
                case "transport":
                    distill_transport = 1
                    break
            }
        }

        if (distill_carbon == 1) {
            ch_distill_carbon = file(params.distill_carbon_sheet).exists() ? file(params.distill_carbon_sheet) : error("Error: If using --distill_topic carbon (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

        } else{
            ch_distill_carbon = default_channel
        }
        if (distill_energy == 1) {
            ch_distill_energy = file(params.distill_energy_sheet).exists() ? file(params.distill_energy_sheet) : error("Error: If using --distill_topic energy (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

        } else{
            ch_distill_energy = default_channel
        }

        if (distill_misc == 1) {
            ch_distill_misc = file(params.distill_misc_sheet).exists() ? file(params.distill_misc_sheet) : error("Error: If using --distill_topic misc (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

        } else{
            ch_distill_misc = default_channel
        }

        if (distill_nitrogen == 1) {
            ch_distill_nitrogen = file(params.distill_nitrogen_sheet).exists() ? file(params.distill_nitrogen_sheet) : error("Error: If using --distill_topic nitrogen (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

        } else{
            ch_distill_nitrogen = default_channel
        }

        if (distill_transport == 1) {
            ch_distill_transport = file(params.distill_transport_sheet).exists() ? file(params.distill_transport_sheet) : error("Error: If using --distill_topic transport (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

        } else{
            ch_distill_transport = default_channel
        }
    } 
    
    if (params.distill_ecosystem != "") {
        distill_eng_sys = 0
        distill_ag = 0

        def distillEcosystemList = params.distill_ecosystem.split()

        // Create a list to store the generated channels
        def ecoSysChannels = []

        distillEcosystemList.each { ecosysItem ->
            if (!['eng_sys', 'ag'].contains(ecosysItem)) {
                error("Invalid distill ecosystem: $ecosysItem. Valid values are eng_sys, ag")
            }

            switch (ecosysItem) {
                case "ag":
                    distill_ag = 1
                    break
                case "eng_sys":
                    distill_eng_sys = 1
                    break
            }
        }

        if (distill_eng_sys == 1) {
            ch_distill_eng_sys = file(params.distill_eng_sys_sheet).exists() ? file(params.distill_eng_sys_sheet) : error("Error: If using --distill_ecosystem eng_sys, you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

        } else{
            ch_distill_eng_sys = default_channel
        }

        if (distill_ag == 1) {
            ch_distill_ag = file(params.distill_ag_sheet).exists() ? file(params.distill_ag_sheet) : error("Error: If using --distill_ecosystem ag, you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")

        } else{
            ch_distill_ag = default_channel
        }
    }
    else{
        ch_distill_ecosys = default_channel
    }
    
    if (params.distill_custom != "") {
        ch_distill_custom = file(params.distill_custom).exists() ? file(params.distill_custom) : error("Error: If using --distill_custom <path/to/TSV>, you must have the preformatted custom distill sheet in the provided file: ${params.distill_custom}.")
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
if( params.call || params.annotate || params.distill_ecosystem || params.distill_topic){
    log.info """
            DRAM2 Nextflow
            ===================================
            fastas       : ${params.input_fasta}
            outdir       : ${params.outdir}
            threads      : ${params.threads}
            call genes   : ${params.call ? 'true' : 'false'}
            annotate     : ${params.annotate ? 'true' : 'false'}
            distill      : ${params.distill ? 'true' : 'false'}
            """
            .stripIndent()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow for pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    /* Rename fasta headers
        Process 1-by-1 
    */
    if( params.annotate || params.call ){
        if( params.rename ) {
            RENAME_FASTA( fastas )
            fasta = RENAME_FASTA.out.renamed_fasta
        } 
        else {
            fasta = fastas
        }
    }

    /* Call genes using prodigal - only if the user did not provide input genes */
    if( params.call && params.input_genes == 0 ) {
        CALL_GENES ( fasta )
        called_genes = CALL_GENES.out.prodigal_fna
        called_proteins = CALL_GENES.out.prodigal_faa
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Annotation
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if( params.annotate ){
        // Not sure if we want these default or optional yet
        /* Run tRNAscan-SE on each fasta to identify tRNAs */
        TRNA_SCAN( fasta )
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
        RRNA_SCAN( fasta )
        ch_rrna_scan = RRNA_SCAN.out.rrna_scan_out
        Channel.empty()
            .mix( ch_rrna_scan )
            .collect()
            .set { ch_collected_rRNAs }
        /* Run RRNA_COLLECT to generate a combined TSV for all fastas */
        RRNA_COLLECT( ch_collected_rRNAs )
        ch_rrna_sheet = RRNA_COLLECT.out.rrna_collected_out
        ch_rrna_combined = RRNA_COLLECT.out.rrna_combined_out


        /* Annotate according to the user-specified databases */
        if( annotate_kegg == 1 ){
            //KEGG_INDEX ( params.kegg_mmseq_loc )
            //MMSEQS2 ( called_genes, params.kegg_mmseq_loc, params.kegg_index )
            //MMSEQS2 ( called_genes, params.kegg_mmseq_loc )
        } else {
            ch_kegg_formatted = []
        }
        if( annotate_kofam == 1 ){
            HMM_SEARCH_KOFAM ( called_proteins, ch_kofam_db )
            ch_kofam_hmms = HMM_SEARCH_KOFAM.out.hmm_search_out

            PARSE_HMM_KOFAM ( ch_kofam_hmms, ch_parse_hmmsearch )
            ch_kofam_parsed = PARSE_HMM_KOFAM.out.parsed_hmm

            KOFAM_HMM_FORMATTER ( ch_kofam_parsed, params.kofam_top_hit, ch_kofam_list, ch_kofam_formatter )
            ch_kofam_formatted = KOFAM_HMM_FORMATTER.out.kofam_formatted_hits
        } else {
            ch_kofam_formatted = []
        }
        if( annotate_dbcan == 1 ){
            
            HMM_SEARCH_DBCAN ( called_proteins, ch_dbcan_db )
            ch_dbcan_hmms = HMM_SEARCH_DBCAN.out.hmm_search_out

            PARSE_HMM_DBCAN ( ch_dbcan_hmms, ch_parse_hmmsearch )
            ch_dbcan_parsed = PARSE_HMM_DBCAN.out.parsed_hmm

            DBCAN_HMM_FORMATTER ( ch_dbcan_parsed, params.dbcan_top_hit, ch_dbcan_fam, ch_dbcan_subfam, ch_dbcan_formatter )
            ch_dbcan_formatted = DBCAN_HMM_FORMATTER.out.dbcan_formatted_hits
            
        } else {
            ch_dbcan_formatted = []
        }
        if (annotate_camper == 1){
        }
        if (annotate_fegenie == 1){
        }
        if (annotate_methyl == 1){
        }
        if (annotate_cant_hyd == 1){
        }
        if (annotate_heme == 1){
        }
        if (annotate_sulfur == 1){
        }
        if (annotate_methyl == 1){

        }
        /*
        if(user provided database){

        }
        */

        /* Combine formatted annotations */
        /* Collect all sample formatted_hits in prep for distill_summary */
        // Need to figure out how to handle when not all channels are here.
        Channel.empty()
            .mix( ch_kofam_formatted )
            .mix( ch_dbcan_formatted )
            .collect()
            .set { collected_formatted_hits }

        /* COMBINE_ANNOTATIONS collects all annotations files across ALL databases */
        COMBINE_ANNOTATIONS( collected_formatted_hits, ch_combine_annot_script )
        ch_combined_annotations = COMBINE_ANNOTATIONS.out.combined_annotations_out

        COUNT_ANNOTATIONS ( ch_combined_annotations, ch_count_annots_script )
        ch_annotation_counts = COUNT_ANNOTATIONS.out.target_id_counts

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
            ch_final_annots = ADD_TAXA.out.annots_taxa_out
        }
        else{
            ch_final_annots = ch_combined_annotations
        }

    }
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Distill
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */   
    if( params.distill_topic != "" || params.distill_ecosystem != "" || params.distill_custom != "" )
    {
        /* Combine the individual user-specified distill sheets into a single channel */
        COMBINE_DISTILL(ch_distill_carbon, ch_distill_energy, ch_distill_misc, ch_distill_nitrogen, ch_distill_transport, ch_distill_ag, ch_distill_eng_sys, ch_distill_custom )
        ch_combined_distill_sheets = COMBINE_DISTILL.out.ch_combined_distill_sheets

        DISTILL_SUMMARY( ch_final_annots, ch_combined_distill_sheets, ch_annotation_counts, ch_distill_summary_script )
        ch_simple_matab_summ = DISTILL_SUMMARY.out.ch_genome_sum_simple

        //Need to add in distill final which make the multi-sheet xlsx:
        // 1) add in functionality to process Bin Quality and Taxonomy (if present on the ch_final_annots channel)
        DISTILL_FINAL( ch_simple_matab_summ, ch_distill_final_script, ch_rrna_sheet, ch_trna_sheet, ch_final_annots )
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
    Input helper-function definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def validOptions = ["--call", "--annotate", "--distill"]


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Define a function to build additional parameters string for DRAM distill sheets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



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