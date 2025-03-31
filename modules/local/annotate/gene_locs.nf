process GENE_LOCS {

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2:89f055454dac3575"
    
    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( genes )
    path( ch_called_genes_loc_script_faa )

    output:
    tuple val( input_fasta ), path( "${input_fasta}_called_genes_table.tsv" ), emit: prodigal_locs_tsv


    script:

    """
    generate_faa_gene_loc_tsv.py ${genes} "${input_fasta}_called_genes_table.tsv"
    """
}
