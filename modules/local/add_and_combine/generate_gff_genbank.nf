process GENERATE_GFF_GENBANK {
    label 'process_tiny'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_biopython:7df21d027f67112e"
    
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
    # export constants for script
    export FASTA_COLUMN="${params.CONSTANTS.FASTA_COLUMN}"

    mkdir -p GFF
    mkdir -p GBK

    generate_gff_genbank.py ${flags.join(' ')} --input_fastas_paths ${all_genes_fna} --database_list ${databases_list} --annotations ${raw_annotations}
    """
}