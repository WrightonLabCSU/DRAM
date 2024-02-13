process MMSEQS_INDEX{

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