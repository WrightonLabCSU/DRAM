process MMSEQS_INDEX{
    label 'process_tiny'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:fc59f737a5e0566a"
    
    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( fasta )

    output:
    tuple val( input_fasta ), path( "*.mmsdb*" ), emit: mmseqs_index_out

    script:
    """
    # Create temporary directory
    mkdir -p mmseqs_out/tmp

    mmseqs createdb ${fasta} ${input_fasta}.mmsdb
    mmseqs createindex ${input_fasta}.mmsdb mmseqs_out/tmp --threads ${params.threads}

    """


}