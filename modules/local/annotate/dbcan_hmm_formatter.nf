process DBCAN_HMM_FORMATTER {
    
    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    
    tag { sample }

    input:
    tuple val( sample ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    val(db_name)
    file(ch_sql_descriptions_db)

    output:
    tuple val( sample ), path ( "${sample}_formatted_dbcan_hits.csv" ), emit: dbcan_formatted_hits
    tuple val(sample), path("${sample}_sql_formatted_${db_name}_hits.csv"), emit: sql_formatted_hits


    script:
    """
    dbcan_hmm_formatter.py --hits_csv ${hits_file} --output "${sample}_formatted_dbcan_hits.csv" --gene_locs "gene_locs.tsv"

    sql_add_descriptions.py --hits_csv "${sample}_formatted_dbcan_hits.csv" --db_name ${db_name} --output "${sample}_sql_formatted_${db_name}_hits.csv" --db_file ${ch_sql_descriptions_db}

    
    """
}
