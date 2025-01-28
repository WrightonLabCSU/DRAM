process GENE_LOCS {

    tag { sample }

    input:
    tuple val( sample ), path( genes )
    path( ch_called_genes_loc_script_faa )

    output:
    tuple val( sample ), path( "${sample}_called_genes_table.tsv" ), emit: prodigal_locs_tsv


    script:

    """


    python ${ch_called_genes_loc_script_faa} ${genes} "${sample}_called_genes_table.tsv"

    """
}
