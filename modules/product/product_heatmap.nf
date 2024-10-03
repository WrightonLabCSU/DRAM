/*
This is purely a placeholder process for future Product development and can be modified freely
*/
process PRODUCT_HEATMAP {

    input:
    path( ch_final_annots, stageAs: "raw-annotations.tsv")
    val(groupby_column)

    output:
    path( "product.html" ), emit: product_html
    path( "product.tsv" ), emit: product_tsv

    script:
    """
    python -m ${params.make_product_pkg} \\
    --annotations ${ch_final_annots} \\
    --groupby-column ${groupby_column} \\

    """
}
