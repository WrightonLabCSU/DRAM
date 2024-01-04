process DBCAN_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file )
    val( top_hit )
    file( ch_dbcan_formatter )

    output:
    tuple val( sample ), path ( "${sample}_formatted_hits.out" ), emit: dbcan_formatted_hits


    output:
    path ( "formatted_dbcan_hits.out" ), emit: formatted_hits

    script:
    """
    python ${ch_dbcan_formatter} --hits_csv ${hits_file} --db_name ${database_name} --db_handler ${db_handler_file} --output "formatted_dbcan_hits.out"
    
    """
}
