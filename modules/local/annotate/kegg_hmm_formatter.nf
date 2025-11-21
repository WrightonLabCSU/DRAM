process KEGG_HMM_FORMATTER {
    label 'process_tiny'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:d2c88b719ab1322c"

    tag { input_fasta }

    input:
    file( hmm_info_path )
    tuple val( input_fasta ), path( hits_file )
    val( top_hit )

    output:
    tuple val( input_fasta ), path ( "${input_fasta}___formatted_kegg_hits.csv" ), emit: formatted_hits

    script:
    """
    kegg_hmm_formatter.py --hits_csv ${hits_file} --hmm_info_path ${hmm_info_path} --top_hit "${top_hit}" --output "${input_fasta}___formatted_kegg_hits.csv"
    """
}
