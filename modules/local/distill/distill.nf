process DISTILL {
    label 'process_medium'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_openpyxl:20d9b4833f4ee53c"

    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )
    path( ch_combined_distill_sheets, stageAs: "combined/*" )
    path( ch_target_id_counts, stageAs: "target_id_counts.tsv")
    path( ch_quast_stats, stageAs: "collected_quast.tsv" )
    path( ch_rrna_sheet, stageAs: "rrna_sheet.tsv" )
    path( ch_combined_rrna, stageAs: "rrna_combined.tsv" )
    path( ch_trna_sheet, stageAs: "trna_sheet.tsv" )
    path( annotations_sqlite3 )

    output:
    path( "distillate.xlsx" ), emit: distillate

    script:
    """
    # export constants for script
    export FASTA_COLUMN="${params.CONSTANTS.FASTA_COLUMN}"

    distill_xlsx.py --target_id_counts ${ch_target_id_counts} --db_name ${annotations_sqlite3} --distill_sheets combined/*.tsv --rrna_file ${ch_rrna_sheet} --combined_rrna_file ${ch_combined_rrna} --trna_file ${ch_trna_sheet} --quast ${ch_quast_stats} --output_file "distillate.xlsx" --threads ${params.threads}

    """
}
