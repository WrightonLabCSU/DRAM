MMSEQS_INDEX{
    tag { sample }

    input:
    tuple val( sample ), path( fasta )

    output:
    tuple val( sample ), path( "${sample}_search_results.tsv" ), emit: mmseqs_index_out

    script:
    """

    mmseqs createdb ${fasta} ${sample}.mmsdb
    
    """


}