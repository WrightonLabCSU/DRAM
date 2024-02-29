process DISTILL {

    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )
    path( ch_combined_distill_sheets )
    path( ch_annotation_counts, stageAs: "target_id_counts.tsv")
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
    
   
    """
}

//python ${ch_distill_xlsx_script}  --db_name annotations.db --distill_sheets ${distill_sheets} --rrna_file ${rrna_sheet} --trna_file ${trna_sheet} --output_file "distillate.xlsx"