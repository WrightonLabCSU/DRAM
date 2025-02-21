process MMSEQS_INDEX{

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    
    tag { sample }

    input:
    tuple val( sample ), path( fasta )

    output:
    tuple val( sample ), path( "*.mmsdb*" ), emit: mmseqs_index_out

    script:
    """

    mmseqs createdb ${fasta} ${sample}.mmsdb

    """


}