process CALL_GENES {

    tag { sample }

    input:
    tuple val( sample ), path( fasta )
    path( ch_called_genes_loc_script )

    output:
    tuple val( sample ), path( "${sample}_called_genes.fna" ), emit: prodigal_fna
    tuple val( sample ), path( "${sample}_called_genes.faa" ), emit: prodigal_faa
    tuple val( sample ), path( "${sample}_called_genes_table.tsv" ), emit: prodigal_locs_tsv
    path( "${sample}_${params.min_contig_len}.fa" ), emit: prodigal_filtered_fasta
    path( "${sample}_called_genes.gff" ), emit: prodigal_gff


    script:

    """

    /opt/bbmap/reformat.sh \\
    in=${fasta} \\
    out="${sample}_${params.min_contig_len}.fa" \\
    minlength=${params.min_contig_len}

    prodigal \\
    -i "${sample}_${params.min_contig_len}.fa" \\
    -o "${sample}_called_genes.gff" \\
    -p ${params.prodigal_mode} \\
    -g ${params.prodigal_trans_table} \\
    -f gff \\
    -d "${sample}_called_genes.fna" \\
    -a "${sample}_called_genes.faa"

    python ${ch_called_genes_loc_script} "${sample}_called_genes.gff" > "${sample}_called_genes_table.tsv"

    """
}
