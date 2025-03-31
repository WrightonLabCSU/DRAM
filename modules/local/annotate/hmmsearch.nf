process HMM_SEARCH {

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2:89f055454dac3575"

    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( fasta )
    val ( e_value )
    path( database_loc )

    output:
    tuple val( input_fasta ), path("${input_fasta}_hmmsearch.out"), emit: hmm_search_out

    script:

    """
    ln -s ${database_loc}/* .

    hmmsearch \\
    -E ${e_value} \\
    --domtblout ${input_fasta}_hmmsearch.out \\
    --cpu 2 \\
    *.hmm \\
    ${fasta}

    """

}
