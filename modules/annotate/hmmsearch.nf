process HMM_SEARCH {

    tag { sample }

    input:
    tuple val( sample ), path( fasta )
    path( database_loc )

    output:
    tuple val( sample ), path("${sample}_hmmsearch.out"), emit: hmm_search_out

    script:

    """

    hmmsearch \\
    --domtblout ${sample}_hmmsearch.out \\
    --cpu 2 \\
    "${database_loc}/*" \\
    ${fasta}

    """

}
