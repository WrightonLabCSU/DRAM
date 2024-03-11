process MMSEQS_SEARCH {

    errorStrategy 'finish'

    tag { sample }

    input:
    tuple val( sample ), path( fasta ), path(gff)

    output:
    tuple val( sample ), path( "${sample}_called_genes.fna" ), emit: prodigal_fna

    script:
    """
    quast.py \\
    ${fasta} \\ 
    -g ${gff} \\
    -o "${sample}_QUAST" \\
    --threads ${params.threads}

    """

}