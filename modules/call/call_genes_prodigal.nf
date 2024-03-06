process CALL_GENES {

    tag { sample }

    input:
    tuple val( sample ), path( fasta )
    path( ch_called_genes_loc_script )

    output:
    tuple val( sample ), path( "${sample}_called_genes.fna" ), emit: prodigal_fna
    tuple val( sample ), path( "${sample}_called_genes.faa" ), emit: prodigal_faa
    tuple val( sample ), path( "${sample}_called_genes.${params.prodigal_format}" ), emit: prodigal_output
    tuple val( sample ), path( "${sample}_called_genes_table.tsv" ), emit: prodigal_locs_tsv


    script:

    """

    prodigal \\
    -i ${fasta} \\
    -o "${sample}_called_genes.${params.prodigal_format}" \\
    -p ${params.prodigal_mode} \\
    -g ${params.prodigal_trans_table} \\
    -f ${params.prodigal_format} \\
    -d "${sample}_called_genes.fna" \\
    -a "${sample}_called_genes.faa"

    python ${ch_called_genes_loc_script} "${sample}_called_genes.${params.prodigal_format}" > "${sample}_called_genes_table.tsv"

    """
}
