process ADD_SQL_DESCRIPTIONS {

    errorStrategy 'finish'

    input:
    tuple val(sample), path(hits_file)
    val(db_name)
    file(ch_sql_descriptions_db)
    file(ch_sql_parser)

    output:
    tuple val(sample), path("${sample}_sql_formatted_${db_name}_hits.csv"), emit: sql_formatted_hits

    script:
    """
    python ${ch_sql_parser} --hits_csv ${hits_file} --db_name ${db_name} --output "${sample}_sql_formatted_${db_name}_hits.csv" --db_file ${ch_sql_descriptions_db}
    """
}
