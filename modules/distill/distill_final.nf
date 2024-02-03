process DISTILL_FINAL {

    input:
    path( metabolism_summary )
    file( ch_distill_final_script )
    file( ch_rrna_sheet, stageAs: "rrna_sheet.tsv" )
    file( combined_rrna, stageAs: "rrna_combined.tsv" )
    file( ch_trna_sheet, stageAs: "trna_sheet.tsv" )
    file( combined_annotations )

    output:
    path( "distillate.xlsx" ), emit: distillate

    script:
    """

    python ${ch_distill_final_script} --input_file ${metabolism_summary} --rrna_file ${ch_rrna_sheet} --combined_rrna ${combined_rrna} --trna_file ${ch_trna_sheet} --combined_annotations ${combined_annotations} --output_file distillate.xlsx

    """
}