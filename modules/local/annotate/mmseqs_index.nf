process MMSEQS_INDEX{

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2:89f055454dac3575"
    
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