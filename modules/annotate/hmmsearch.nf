process HMM_SEARCH {

    tag { sample }

    input:
    tuple val( sample ), path( fasta )
    val ( e_value )
    path( database_loc )

    output:
    tuple val( sample ), path("${sample}_hmmsearch.out"), emit: hmm_search_out

    script:

    """
    ln -s ${database_loc}/* .

    hmmsearch \\
    -E ${e_value} \\
    --domtblout ${sample}_hmmsearch.out \\
    --cpu 2 \\
    *.hmm \\
    ${fasta}

    """

}
