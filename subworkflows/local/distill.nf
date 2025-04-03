//
// Subworkflow with functionality specific to the WrightonLabCSU/dram pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DISTILL as DISTILL_SCRIPT                         } from "${projectDir}/modules/local/distill/distill.nf"
include { COMBINE_DISTILL                                   } from "${projectDir}/modules/local/distill/combine_distill.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO DISTILL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DISTILL {
    take:
    ch_distill_carbon
    ch_distill_energy
    ch_distill_misc
    ch_distill_nitrogen
    ch_distill_transport
    ch_distill_ag
    ch_distill_eng_sys
    ch_distill_camper
    ch_distill_custom_collected
    ch_final_annots
    ch_annotation_counts
    ch_rrna_sheet
    ch_quast_stats
    ch_rrna_combined
    ch_trna_sheet
    ch_annotations_sqlite3

    main:

    // Combine the individual user-specified distill sheets into a single channel
    COMBINE_DISTILL(ch_distill_carbon, ch_distill_energy, ch_distill_misc, ch_distill_nitrogen, ch_distill_transport, ch_distill_ag, ch_distill_eng_sys, ch_distill_camper, ch_distill_custom_collected )
    ch_combined_distill_sheets = COMBINE_DISTILL.out.ch_combined_distill_sheets

    // Generate multi-sheet XLSX document containing annotations included in user-specified distillate speadsheets

    DISTILL_SCRIPT( ch_final_annots, ch_combined_distill_sheets, ch_annotation_counts, ch_quast_stats, ch_rrna_sheet, ch_rrna_combined, ch_trna_sheet, ch_annotations_sqlite3 )
    ch_distillate = DISTILL_SCRIPT.out.distillate


}