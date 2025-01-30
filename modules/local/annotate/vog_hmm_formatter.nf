process VOG_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    val( top_hit )
    val( db_name )
    file(ch_sql_descriptions_db)

    output:
    tuple val( sample ), path ( "${sample}_formatted_vog_hits.out" ), emit: vog_formatted_hits
    tuple val(sample), path("${sample}_sql_formatted_${db_name}_hits.csv"), emit: sql_formatted_hits

    script:
    """
    vog_hmm_formatter.py --hits_csv ${hits_file} --gene_locs ${prodigal_locs_tsv} --output "${sample}_formatted_vog_hits.out"

    sql_add_descriptions.py --hits_csv "${sample}_formatted_vog_hits.out" --db_name ${db_name} --output "${sample}_sql_formatted_${db_name}_hits.csv" --db_file ${ch_sql_descriptions_db}
    """
}
