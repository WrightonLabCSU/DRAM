process BOWTIE_DATABASE {

    input:
    path(clustered_fasta)
    path(clustered_files)
    
    output:
    path("*gene_DB*"), emit: gene_db


    script:

    """
    source /opt/miniconda/bin/activate support

    bowtie2-build \\
    ${clustered_fasta} \\
    gene_DB \\
    --large-index \\
    --threads ${params.threads}

    conda deactivate
    
    """
}  