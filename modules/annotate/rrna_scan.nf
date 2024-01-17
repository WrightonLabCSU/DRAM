process RRNA_SCAN {

    errorStrategy 'finish'

    tag { sample }

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}__processed_rrnas.tsv"), emit: trna_scan_out, optinal: true

    script:

    """
    barrnap \\
    --cpu ${params.threads} \\
    ${fasta}



    """

}
