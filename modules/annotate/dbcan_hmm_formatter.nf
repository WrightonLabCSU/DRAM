process DBCAN_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    val( top_hit )
    val(db_name)
    file(ch_dbcan_formatter)
    file(ch_sql_parser)
    file(ch_sql_descriptions_db)

    output:
    tuple val( sample ), path ( "${sample}_formatted_dbcan_hits.csv" ), emit: dbcan_formatted_hits
    tuple val(sample), path("${sample}_sql_formatted_${db_name}_hits.csv"), emit: sql_formatted_hits


    script:
    """
    python ${ch_dbcan_formatter} --hits_csv ${hits_file} --output "${sample}_formatted_dbcan_hits.csv" --gene_locs "gene_locs.tsv"

    python ${ch_sql_parser} --hits_csv "${sample}_formatted_dbcan_hits.csv" --db_name ${db_name} --output "${sample}_sql_formatted_${db_name}_hits.csv" --db_file ${ch_sql_descriptions_db}

    
    """
}

