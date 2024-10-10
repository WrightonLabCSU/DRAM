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
    // Correct usage: Directly use params without ${...} for Groovy code
    def gff_flag = params.generate_gff ? "--gff" : ""
    def gbk_flag = params.generate_gbk ? "--gbk" : ""
    def flags = []
    if(params.generate_gff) flags += "--gff"
    if(params.generate_gbk) flags += "--gbk"

    // Make sure to properly pass the `flags` variable to the Python script
    """
    mkdir -p GFF
    mkdir -p GBK

    python ${ch_generate_gff_gbk} ${flags.join(' ')} --samples_paths ${all_genes_fna} --database_list ${databases_list} --annotations ${raw_annotations}
    """
}