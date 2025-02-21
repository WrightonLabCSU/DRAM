/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_dram_pipeline'
include { getFastaChannel        } from '../subworkflows/local/utils_pipeline_setup.nf'

// Pipeline steps
include { RENAME_FASTA           } from "${projectDir}/modules/local/rename/rename_fasta.nf"
include { CALL                   } from "${projectDir}/subworkflows/local/call.nf"
include { COLLECT_RNA            } from "${projectDir}/subworkflows/local/collect_rna.nf"
include { MERGE                  } from "${projectDir}/subworkflows/local/merge.nf"
include { ANNOTATE               } from "${projectDir}/subworkflows/local/annotate.nf"
include { ADD_AND_COUNT        } from "${projectDir}/subworkflows/local/add_and_count.nf"
include { DISTILL                } from "${projectDir}/subworkflows/local/distill.nf"


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

    default_channel = Channel.fromPath(params.distill_dummy_sheet)
    ch_fasta = getFastaChannel(params.input_fasta, params.fasta_fmt)
    distill_flag = (params.distill_topic != "" || params.distill_ecosystem != "" || params.distill_custom != "")


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

    if (distill_flag) {
        if (params.distill_topic != "") {
            def validTopics = ['default', 'carbon', 'energy', 'misc', 'nitrogen', 'transport', 'camper']
            def topics = params.distill_topic.split(',')

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
            def distillEcosystemList = params.distill_ecosystem.split(',')

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
    // Pipeline steps
    //

    

    if( params.rename ) {
        RENAME_FASTA( ch_fasta )
        ch_fasta = RENAME_FASTA.out.renamed_fasta
    }
    if (params.merge_annotations){
        MERGE()
    }else {
        ch_quast_stats = default_channel
        ch_gene_locs = default_channel
        ch_called_proteins = default_channel
        ch_collected_fna = default_channel

        if (params.call){
            CALL( ch_fasta )
            ch_quast_stats = CALL.out.ch_quast_stats
            ch_gene_locs = CALL.out.ch_gene_locs
            ch_called_proteins = CALL.out.ch_called_proteins
            ch_collected_fna = CALL.out.ch_collected_fna

        } else{
            
        }

        if (params.annotate){
            ANNOTATE( ch_fasta, ch_gene_locs, ch_called_proteins, default_channel )
            
        }

        if (params.call || distill_flag){
            COLLECT_RNA( ch_fasta )
        }
        if (params.annotate || distill_flag){
            if (params.annotate){
                ch_combined_annotations = ANNOTATE.out.ch_combined_annotations
            } else {
                ch_combined_annotations = Channel
                    .fromPath(params.annotations, checkIfExists: true)
                    .ifEmpty { exit 1, "If you specify --distill_<topic|ecosystem|custom> without --annotate, you must provide an annotations TSV file (--annotations <path>) with approprite formatting. Cannot find any called gene files matching: ${params.annotations}\nNB: Path needs to follow pattern: path/to/directory/" }
            }
            ADD_AND_COUNT( ch_combined_annotations )

            
            if( distill_flag ){

                DISTILL(
                    ch_distill_carbon,
                    ch_distill_energy,
                    ch_distill_misc,
                    ch_distill_nitrogen,
                    ch_distill_transport,
                    ch_distill_ag,
                    ch_distill_eng_sys,
                    ch_distill_camper,
                    ch_distill_custom_collected,
                    ADD_AND_COUNT.out.ch_final_annots,
                    ADD_AND_COUNT.out.ch_annotation_counts,
                    COLLECT_RNA.out.ch_rrna_sheet,
                    ch_quast_stats,
                    COLLECT_RNA.out.ch_rrna_combined,
                    COLLECT_RNA.out.ch_trna_sheet,
                    ADD_AND_COUNT.out.ch_annotations_sqlite3
                )
            }
        }
    }

    // //
    // // MODULE: MultiQC
    // //
    // ch_multiqc_config        = Channel.fromPath(
    //     "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // ch_multiqc_custom_config = params.multiqc_config ?
    //     Channel.fromPath(params.multiqc_config, checkIfExists: true) :
    //     Channel.empty()
    // ch_multiqc_logo          = params.multiqc_logo ?
    //     Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
    //     Channel.empty()

    // summary_params      = paramsSummaryMap(
    //     workflow, parameters_schema: "nextflow_schema.json")
    // ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
    //     file(params.multiqc_methods_description, checkIfExists: true) :
    //     file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    // ch_methods_description                = Channel.value(
    //     methodsDescriptionText(ch_multiqc_custom_methods_description))

    // ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_methods_description.collectFile(
    //         name: 'methods_description_mqc.yaml',
    //         sort: true
    //     )
    // )

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList(),
    //     [],
    //     []
    // )

    emit:
    // multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
