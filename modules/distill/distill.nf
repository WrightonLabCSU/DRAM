process DISTILL {

    input:
    file( combined_annotations, stageAs: "raw-annotations.tsv" )
    val( database_names_list )
    path( ch_combined_distill_sheets )
    path( ch_rrna_sheet, stageAs: "rrna_sheet.tsv" )
    path( combined_rrna, stageAs: "rrna_combined.tsv" )
    path( ch_trna_sheet, stageAs: "trna_sheet.tsv" )
    file( ch_distill_xlsx_script )
    file( ch_distill_sql_script )

    output:
    path( "distillate.xlsx" ), emit: distillate

    script:
    """

    python ${ch_distill_sql_script} --combined_annotations ${combined_annotations} --db_name "annotations.db" --db_list ${database_names_list}
    
   
    """
}

//python ${ch_distill_xlsx_script}  --db_name annotations.db --distill_sheets ${distill_sheets} --rrna_file ${rrna_sheet} --trna_file ${trna_sheet} --output_file "distillate.xlsx"