/*
This is purely a placeholder process for future Product development and can be modified freely
*/
process PRODUCT_HEATMAP {

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"

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
