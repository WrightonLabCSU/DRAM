//
// Subworkflow with functionality specific to the WrightonLabCSU/dram pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TRNA_SCAN                                     } from "${projectDir}/modules/local/collect_rna/trna_scan.nf"
include { RRNA_SCAN                                     } from "${projectDir}/modules/local/collect_rna/rrna_scan.nf"
include { TRNA_COLLECT                                  } from "${projectDir}/modules/local/collect_rna/trna_collect.nf"
include { RRNA_COLLECT                                  } from "${projectDir}/modules/local/collect_rna/rrna_collect.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO COLLECT RNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COLLECT_RNA {
    take:
    ch_fasta  // channel: [ val(sample name), path(fasta) ]
    default_channel // channel: dummy sheet

    main:


    // If we didn't run call 
    if (!params.call) {
        // the user provided rrnas or trnas 
        if( params.trnas != "" ){
            Channel.fromPath("${params.trnas}/*.tsv", checkIfExists: true)
                .ifEmpty { exit 1, "If you specify --distill_<topic|ecosystem|custom> without --call, you must provide individual rRNA files generated with tRNAscan-SE. Cannot find any files at: ${params.trnas}\nNB: Path needs to follow pattern: path/to/directory" }
                .collect()
                .set { ch_collected_tRNAs }
        }
        if( params.rrnas != "" ){
            Channel.fromPath("${params.rrnas}/*.tsv", checkIfExists: true)
                .ifEmpty { exit 1, "If you specify --distill_<topic|ecosystem|custom> without --call, you must provide individual rRNA files generated with barrnap. Cannot find any files at: ${params.rrnas}\nNB: Path needs to follow pattern: path/to/directory" }
                .collect()
                .set { ch_collected_rRNAs }
        }
    } else { // If we did run call then we need to generate the rrnas and trnas from the fastas
        // Run tRNAscan-SE on each fasta to identify tRNAs
        TRNA_SCAN( ch_fasta )
        ch_trna_scan = TRNA_SCAN.out.trna_scan_out
        // Collect all sample formatted tRNA files
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

    }

    // Create sheet for trnas from the collected tRNAs or provided tRNAs
    if ((params.call) || (params.trnas != "")) {
        // Run TRNA_COLLECT to generate a combined TSV for all fastas
        TRNA_COLLECT( ch_collected_tRNAs )
        ch_trna_sheet = TRNA_COLLECT.out.trna_collected_out
    } else{
        ch_trna_sheet = default_channel
    }

    // Create sheet for rrnas from the collected rRNAs or provided rRNAs
    if ((params.call) || (params.rrnas != "")) {
        // Run RRNA_COLLECT to generate a combined TSV for all fastas
        RRNA_COLLECT( ch_collected_rRNAs )
        ch_rrna_sheet = RRNA_COLLECT.out.rrna_collected_out
        ch_rrna_combined = RRNA_COLLECT.out.rrna_combined_out
    } else {
        ch_rrna_sheet = default_channel
        ch_rrna_combined = default_channel
    }


    emit:
    ch_rrna_sheet
    ch_rrna_combined
    ch_trna_sheet

}