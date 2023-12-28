process CALL_GENES {

    tag { sample }

    input:
    tuple val( sample ), path( fasta )

    output:
    tuple val( sample ), path( "${sample}_called_genes.fna" ), emit: prodigal_fna
    tuple val( sample ), path( "${sample}_called_genes.faa" ), emit: prodigal_faa
    tuple val( sample ), path( "${sample}_progigal_out.${params.prodigal_format}" ), emit: prodigal_output


    script:

    """

    prodigal \\
    -i ${fasta} \\
    -o "${sample}_progigal_out.${params.prodigal_format}" \\
    -p ${params.prodigal_mode} \\
    -g ${params.prodigal_trans_table} \\
    -f ${params.prodigal_format} \\
    -d "${sample}_called_genes.fna" \\
    -a "${sample}_called_genes.faa"

    """
}
