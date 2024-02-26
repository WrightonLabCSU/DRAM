process DBCAN_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file )
    val( top_hit )
    file(ch_sql_descriptions_db)
    file(ch_dbcan_formatter)

    output:
    tuple val( sample ), path ( "${sample}_formatted_dbcan_hits.out" ), emit: dbcan_formatted_hits


    script:
    """
    python ${ch_dbcan_formatter} --hits_csv ${hits_file} --output "${sample}_formatted_dbcan_hits.out" --db_file ${ch_sql_descriptions_db}
    
    """
}

