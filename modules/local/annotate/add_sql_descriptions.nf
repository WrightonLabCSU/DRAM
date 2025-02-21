process ADD_SQL_DESCRIPTIONS {

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample), path(hits_file)
    val(db_name)
    file(ch_sql_descriptions_db)

    output:
    tuple val(sample), path("${sample}_sql_formatted_${db_name}_hits.csv"), emit: sql_formatted_hits

    script:
    """
    sql_add_descriptions.py --hits_csv ${hits_file} --db_name ${db_name} --output "${sample}_sql_formatted_${db_name}_hits.csv" --db_file ${ch_sql_descriptions_db}
    """
}
