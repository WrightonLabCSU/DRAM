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
    ch_final_annots
    ch_rrna_combined
    ch_trna_combined


    main:

    // Generate multi-sheet XLSX document containing annotations included in user-specified distillate speadsheets

    DISTILL_SCRIPT( ch_final_annots, ch_rrna_combined, ch_trna_combined )
    ch_distillate = DISTILL_SCRIPT.out.distillate


}