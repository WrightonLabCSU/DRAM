process MMSEQS_INDEX{

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    
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