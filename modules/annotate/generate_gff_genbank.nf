process GENERATE_GFF_GENBANK {

    input:
    val( all_genes_fna )
    val( databases_list )
    path( raw_annotations, stageAs: "raw-annotations.tsv" )
    path( ch_generate_gff_gbk )

    output:
    path( "GFF/" ), emit: output_gff, optional: true
    path( "GBK/" ), emit: output_gbk, optional: true


    script:

    def gff_flag = ${params.generate_gff} ? "--gff" : ""
    def gbk_flag = ${params.generate_gbk} ? "--gbk" : ""

    """
    mkdir -p GFF
    mkdir -p GBK

    python ${ch_generate_gff_gbk} ${flags.join(' ')} --database_list ${database_list} --annotations raw-annotations.tsv
    """
}