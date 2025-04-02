process PRODUCT_HEATMAP {
    label 'process_low'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_dram-viz:8399ca59227b4fd9"

    input:
    path( ch_final_annots, stageAs: "raw-annotations.tsv")
    val(groupby_column)

    output:
    path( "product.html" ), emit: product_html
    path( "product.tsv" ), emit: product_tsv

    script:
    """
    dram_viz --annotations ${ch_final_annots} --groupby-column ${groupby_column}

    """
}
