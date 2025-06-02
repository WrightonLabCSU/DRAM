process DBCAN_HMM_FORMATTER {
    label 'process_tiny'
    
    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:7b4f27307e83be0e"
    
    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    val(db_name)
    file(ch_sql_descriptions_db)

    output:
    tuple val( input_fasta ), path ( "${input_fasta}_formatted_dbcan_hits.csv" ), emit: dbcan_formatted_hits
    tuple val(input_fasta), path("${input_fasta}_sql_formatted_${db_name}_hits.csv"), emit: sql_formatted_hits


    script:
    """
    dbcan_hmm_formatter.py --hits_csv ${hits_file} --output "${input_fasta}_formatted_dbcan_hits.csv" --gene_locs "gene_locs.tsv"

    sql_add_descriptions.py --hits_csv "${input_fasta}_formatted_dbcan_hits.csv" --db_name ${db_name} --output "${input_fasta}_sql_formatted_${db_name}_hits.csv" --db_file ${ch_sql_descriptions_db}

    
    """
}
