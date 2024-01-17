process RRNA_SCAN {

    errorStrategy 'finish'

    tag { sample }

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}__processed_rrnas.tsv"), emit: rrna_scan_out, optional: true

    script:

    """
    barrnap \\
    --threads ${params.threads} \\
    ${fasta} \\
    > ${sample}_rrnas.tsv



    """

}
