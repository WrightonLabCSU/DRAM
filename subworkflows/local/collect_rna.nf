//
// Subworkflow with functionality specific to the WrightonLabCSU/dram pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TRNA_SCAN                                     } from "../../modules/local/collect_rna/trna_scan.nf"
include { RRNA_SCAN                                     } from "../../modules/local/collect_rna/rrna_scan.nf"
include { TRNA_COLLECT                                  } from "../../modules/local/collect_rna/trna_collect.nf"
include { RRNA_COLLECT                                  } from "../../modules/local/collect_rna/rrna_collect.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO COLLECT RNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COLLECT_RNA {
    take:
    ch_fasta  // channel: [ val(input_fasta name), path(fasta) ]
    default_sheet // Path to dummy sheet
    call          // boolean: whether gene calling flag is set

    main:

    run_rrna_collect = false
    run_trna_collect = false

    // If we didn't run call 
    if (!call) {
        if (params.rrnas) {
            Channel.fromPath("${params.rrnas}/*.tsv", checkIfExists: true)
                .ifEmpty { exit 1, "If you specify --distill_<topic|ecosystem|custom> without --call, you must provide individual rRNA files generated with barrnap. Cannot find any files at: ${params.rrnas}\nNB: Path needs to follow pattern: path/to/directory" }
                .collect()
                .set { ch_collected_rRNAs }
                run_rrna_collect = true

        }
        if (params.trnas) {
            // the user provided rrnas or trnas 
            Channel.fromPath("${params.trnas}/*.tsv", checkIfExists: true)
                .ifEmpty { exit 1, "If you specify --distill_<topic|ecosystem|custom> without --call, you must provide individual tRNA files generated with tRNAscan-SE. Cannot find any files at: ${params.trnas}\nNB: Path needs to follow pattern: path/to/directory" }
                .collect()
                .set { ch_collected_tRNAs }
                run_trna_collect = true
        }
    } else { // If we did run call then we need to generate the rrnas and trnas from the fastas
        // Run tRNAscan-SE on each fasta to identify tRNAs
        TRNA_SCAN( ch_fasta )
        ch_trna_scan = TRNA_SCAN.out.trna_scan_out
        // Collect all input_fasta formatted tRNA files
        Channel.empty()
            .mix( ch_trna_scan )
            .collect()
            .set { ch_collected_tRNAs }
        // Run barrnap on each fasta to identify rRNAs
        RRNA_SCAN( ch_fasta )
        ch_rrna_scan = RRNA_SCAN.out.rrna_scan_out
        Channel.empty()
            .mix( ch_rrna_scan )
            .collect()
            .set { ch_collected_rRNAs }
        run_rrna_collect = true
        run_trna_collect = true
    }

    if (run_rrna_collect) {
        // Create sheet for rrnas from the collected rRNAs or provided rRNAs
        // Run RRNA_COLLECT to generate a combined TSV for all fastas
        RRNA_COLLECT( ch_collected_rRNAs )
        ch_rrna_collected = RRNA_COLLECT.out.rrna_collected_out
        ch_rrna_combined = RRNA_COLLECT.out.rrna_combined_out
    } else {
        ch_rrna_collected = default_sheet
        ch_rrna_combined = default_sheet
    }
    if (run_trna_collect) {
        // Create sheet for trnas from the collected tRNAs or provided tRNAs
        // Run TRNA_COLLECT to generate a combined TSV for all fastas
        TRNA_COLLECT( ch_collected_tRNAs )
        ch_trna_collected = TRNA_COLLECT.out.trna_collected_out
        ch_trna_combined = TRNA_COLLECT.out.trna_combined_out
    } else {
        ch_trna_collected = default_sheet
        ch_trna_combined = default_sheet
    }


    emit:
    ch_rrna_collected
    ch_rrna_combined
    ch_trna_collected
    ch_trna_combined

}