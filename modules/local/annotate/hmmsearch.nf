process HMM_SEARCH {
    label 'process_small'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:d2c88b719ab1322c"

    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( fasta ), path( prodigal_locs_tsv )
    val ( e_value )
    path( database_loc )
    path( hmm_info_path )
    val (ec_from_info )
    val (db_name)

    output:
    tuple val( input_fasta ), path ( "${input_fasta}___formatted_${db_name}_hits.csv" ), emit: formatted_hits, optional: true

    script:
    def args = task.ext.args ?: ""
    def ec_flag = ec_from_info ? "--ec_from_info" : ""

    """
    hmmsearch \\
    -E ${e_value} \\
    --domtblout ${input_fasta}_hmmsearch.out \\
    --cpu ${task.cpus} \\
    ${database_loc}/*.hmm \\
    ${fasta} > /dev/null

    hmm_parser.py \\
        --hmm_domtbl ${input_fasta}_hmmsearch.out \\
        --hmm_info_path ${hmm_info_path} \\
        ${ec_flag} \\
        --gene_locs ${prodigal_locs_tsv} \\
        --db_name ${db_name} \\
        --output "${input_fasta}___formatted_${db_name}_hits.csv"
    """
}
