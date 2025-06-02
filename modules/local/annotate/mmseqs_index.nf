process MMSEQS_INDEX{
    label 'process_tiny'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:7b4f27307e83be0e"
    
    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( fasta )

    output:
    tuple val( input_fasta ), path( "*.mmsdb*" ), emit: mmseqs_index_out

    script:
    """

    mmseqs createdb ${fasta} ${input_fasta}.mmsdb

    """


}