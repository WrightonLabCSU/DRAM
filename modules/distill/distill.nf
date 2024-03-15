process DISTILL {

    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )
    path( ch_combined_distill_sheets )
    path( ch_target_id_counts, stageAs: "target_id_counts.tsv")
    path( ch_quast_stats, stageAs: "collected_quast.tsv" )
    path( ch_rrna_sheet, stageAs: "rrna_sheet.tsv" )
    path( ch_combined_rrna, stageAs: "rrna_combined.tsv" )
    path( ch_trna_sheet, stageAs: "trna_sheet.tsv" )
    path( ch_distill_xlsx_script )
    path( ch_distill_sql_script )

    output:
    path( "distillate.xlsx" ), emit: distillate

    script:
    """

    python ${ch_distill_sql_script} --combined_annotations ${ch_combined_annotations} --db_name "annotations.db" 
    python ${ch_distill_xlsx_script} --target_id_counts ${ch_target_id_counts} --db_name "annotations.db" --distill_sheets ${ch_combined_distill_sheets} --rrna_file ${ch_rrna_sheet} --combined_rrna_file ${ch_combined_rrna} --trna_file ${ch_trna_sheet} --quast ${ch_quast_stats} --output_file "distillate.xlsx" --threads ${params.threads}

    """
}
