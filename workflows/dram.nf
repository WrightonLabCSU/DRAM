/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_dram_pipeline'
include { getDBFlag               } from '../subworkflows/local/utils_pipeline_setup.nf'

// Pipeline steps
include { ADJECTIVES             } from "../modules/local/adjectives/adjectives.nf"
include { PRODUCT_HEATMAP        } from "../modules/local/product/product_heatmap.nf"
include { CAT_KEGG_PEP           } from "../modules/local/database/cat_kegg_pep.nf"
include { FORMAT_KEGG_DB         } from "../modules/local/database/format_kegg_db.nf"
include { MERGE                  } from "../subworkflows/local/merge.nf"
include { ANNOTATE               } from "../subworkflows/local/annotate.nf"
include { ADD_ANNOTATIONS        } from "../modules/local/add_and_combine/add_annotations.nf"
include { ADD_AND_COUNT          } from "../subworkflows/local/add_and_count.nf"
include { SUMMARIZE              } from "../modules/local/distill/distill.nf"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DRAM {

    main:

    //
    // Setup
    //

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_fasta = Channel.empty()

    default_sheet = file(params.distill_dummy_sheet)
    distill_flag = (params.summarize || params.distill_topic != "" || params.distill_ecosystem != "" || params.distill_custom != "" || params.sum_ecos != "")

    // if annotate with raw fasta but no call, we can infer we need to call genes, so set call to true
    // Also, if call is specified, set call to true
    call = false
    if ((params.annotate && params.input_fasta != "") || params.call) {
        call = true
    }
    visualize = false
    if (params.product || params.visualize) {
        visualize = true
    }
    traits = false
    if (params.adjectives || params.traits) {
        traits = true
    }


    if (params.rename || call) {
        ch_fasta = Channel
            .fromPath(file(params.input_fasta) / params.fasta_fmt, checkIfExists: true)
                .ifEmpty { exit 1, "Cannot find any fasta files matching: ${params.input_fasta}\nNB: Path needs to follow pattern: path/to/directory/" }

        ch_fasta = ch_fasta.map {
            fasta_name = it.getBaseName()
            tuple(fasta_name, it)
        }
    }

    use_kegg = params.use_kegg
    use_kofam = params.use_kofam
    use_dbcan = params.use_dbcan
    use_camper = params.use_camper
    use_fegenie = params.use_fegenie
    use_methyl = params.use_methyl
    use_canthyd = params.use_canthyd
    use_sulfur = params.use_sulfur
    use_pfam = params.use_pfam
    use_merops = params.use_merops
    use_uniref = params.use_uniref
    use_vog = params.use_vog

    if (params.anno_dbs != "") {
        anno_dbs = params.anno_dbs.tokenize(',').collect { it.trim().toLowerCase() }
        value_for_all = 'all'
        use_kegg = getDBFlag(anno_dbs, 'kegg', value_for_all)
        use_kofam = getDBFlag(anno_dbs, 'kofam', value_for_all)
        use_dbcan = getDBFlag(anno_dbs, 'dbcan', value_for_all)
        use_camper = getDBFlag(anno_dbs, 'camper', value_for_all)
        use_fegenie = getDBFlag(anno_dbs, 'fegenie', value_for_all)
        use_methyl = getDBFlag(anno_dbs, 'methyl', value_for_all)
        use_canthyd = getDBFlag(anno_dbs, 'canthyd', value_for_all)
        use_sulfur = getDBFlag(anno_dbs, 'sulfur', value_for_all)
        // use_pfam = getDBFlag(anno_dbs, 'pfam', value_for_all)
        // PFAM database is currently disabled in this pipeline due to a bug in the DRAM2 implementation with the PFAM database. It will be re-enabled in a future release.
        use_merops = getDBFlag(anno_dbs, 'merops', value_for_all)
        use_uniref = getDBFlag(anno_dbs, 'uniref', value_for_all)
        use_vog = getDBFlag(anno_dbs, 'vog', value_for_all)
    }



    distill_ecosystem = params.sum_ecos
    if (distill_ecosystem == "") {
        distill_ecosystem = params.distill_ecosystem
    }
    if (params.sum_custom != "") {
        distill_custom = file(params.sum_custom)
    }
    else if (params.distill_custom != "") {
        distill_custom = file(params.distill_custom)
    }
    else {
        distill_custom = ""
    }
    distill_topic = params.distill_topic
    if (distill_flag) {
        if (distill_topic != "") {
            def validTopics = ['default', 'carbon', 'energy', 'misc', 'nitrogen', 'transport', 'camper', 'none']
            def topics = distill_topic.split(',')

            topics.each { topic ->
                if (!validTopics.contains(topic)) {
                    error("Invalid distill topic: $topic. Valid values are ${validTopics.join(',')}. If you included those, try comma separating them without spaces.")
                }

                // Handle the 'default' case by setting all default topics to "1"
            }
        }

        if (distill_ecosystem != "") {
            def validEcos = ['eng_sys', 'ag']
            def distillEcosystemList = distill_ecosystem.split(',')

            distillEcosystemList.each { ecosysItem ->
                if (!validEcos.contains(ecosysItem)) {
                    error("Invalid distill ecosystem: $ecosysItem. Valid values are ${validEcos.join(',')}. If you included those, try comma separating them without spaces.")
                }

            }
        }

        if (distill_custom != "") {
            // Verify the directory exists
            def custom_distill_dir = file(distill_custom)
            if (!custom_distill_dir.exists()) {
                error "Error: The specified custom_distill sheet (--sum_custom or --distill_custom) does not exist: ${distill_custom}"
            }      
        }

        if (params.summarize){
            if (distill_topic == "") {
                distill_topic = "default"
            }
        }

        if (distill_topic == "none") {
            distill_topic = ""
        }

        if (!use_kegg && !use_kofam && !use_dbcan && !use_merops) {
            if (!params.annotations) {
                error("Error: If you are using summarize, you must also use --use_kegg, --use_kofam, --use_dbcan, or --use_merops.")
            }
        }
    }



    //
    // Collate and save software versions
    //

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'dram_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // Single step commands
    //

    if (params.format_kegg){
        if ( params.kegg_pep_root_dir ) {
            CAT_KEGG_PEP( file(params.kegg_pep_root_dir) )
            kegg_pep_f = CAT_KEGG_PEP.out.kegg_pep
        }
        else {
            kegg_pep_f = file(params.kegg_pep_loc)
        }

        if ( !kegg_pep_f.exists() ) {
            error("Error: when using running format_kegg, the kegg_pep_loc file must exists or the user must provide the KEGG pep file root directories. KEGG pep file not found at ${params.kegg_pep_loc} nor kegg pep root directory at ${params.kegg_pep_root_dir}.")
        }

        gene_ko_link_f = params.gene_ko_link_loc && file(params.gene_ko_link_loc).exists() ? file(params.gene_ko_link_loc) : default_sheet
        kegg_download_date = params.kegg_download_date ? params.kegg_download_date : "''"
        skip_gene_ko_link = params.skip_gene_ko_link ? 1 : 0
        FORMAT_KEGG_DB( kegg_pep_f, gene_ko_link_f, kegg_download_date, skip_gene_ko_link )

    } else if (params.merge_annotations){
        MERGE()
    } else {


        //
        // Pipeline steps
        //

        ANNOTATE (
            ch_fasta,
            default_sheet,
            call,
            use_kegg,
            use_kofam,
            use_dbcan,
            use_camper,
            use_fegenie,
            use_methyl,
            use_canthyd,
            use_sulfur,
            use_pfam,
            use_merops,
            use_uniref,
            use_vog
        )

        if (distill_flag) {
            if (params.annotate){ // If the user has specified --annotate, us the outputted annotations
                ch_final_annots = ANNOTATE.out.ch_combined_annotations
                if( params.add_annotations ){
                    ch_add_annots = file(params.add_annotations)
                    ADD_ANNOTATIONS( ANNOTATE.out.ch_combined_annotations, ch_add_annots )
                    ch_final_annots = ADD_ANNOTATIONS.out.combined_annots_out
                }
            } else {  // If the user has not specified --annotate, use the provided annotations
                ch_final_annots = Channel
                    .fromPath(params.annotations, checkIfExists: true)
                    .ifEmpty { exit 1, "If you specify --distill_<topic|ecosystem|custom> without --annotate, you must provide an annotations TSV file (--annotations <path>) with approprite formatting. Cannot find any called gene files matching: ${params.annotations}\nNB: Path needs to follow pattern: path/to/directory/" }
            }

            // ADD_AND_COUNT( ch_combined_annotations, call )  // Do any adding of additional annotations and count the annotations
            // ch_final_annots = ADD_AND_COUNT.out.ch_final_annots


                SUMMARIZE(
                    ch_final_annots,
                    ANNOTATE.out.ch_rrna_collected,
                    ANNOTATE.out.ch_trna_collected,
                    ANNOTATE.out.ch_quast_stats,
                    distill_topic,
                    distill_ecosystem,
                    distill_custom
                )
            

        }
        else if (params.annotations) {
            ch_final_annots = Channel
                .fromPath(params.annotations, checkIfExists: true)
                .ifEmpty { exit 1, "Parameter annotations problem: Cannot find any called gene files matching: ${params.annotations}\nNB: Path needs to follow pattern: path/to/directory/" }
        }

        if (visualize) {  // If the user has specified --product after annotate or distill, generate the product heatmap
            if (!ch_final_annots) {
                error("Error: If you specify --product, you must also specify --annotate or --distill_<topic|ecosystem|custom> to generate the product heatmap or provide an annotations TSV file (--annotations <path>).")
            }
            PRODUCT_HEATMAP( ch_final_annots, params.groupby_column )
        }
        //
        // ADJECTIVES
        //

        if( traits ){
            if (!ch_final_annots) {
                error("Error: If you specify --product, you must also specify --annotate or --distill_<topic|ecosystem|custom> to generate the product heatmap or provide an annotations TSV file (--annotations <path>).")
            }
            ADJECTIVES( ch_final_annots )

        }



    }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
