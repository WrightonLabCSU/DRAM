process GENE_LOCS {

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    
    tag { sample }

    input:
    tuple val( sample ), path( genes )
    path( ch_called_genes_loc_script_faa )

    output:
    tuple val( sample ), path( "${sample}_called_genes_table.tsv" ), emit: prodigal_locs_tsv


    script:

    """
    generate_faa_gene_loc_tsv.py ${genes} "${sample}_called_genes_table.tsv"
    """
}
