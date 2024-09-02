/*
This is purely a placeholder process for future Product development and can be modified freely
*/
process PRODUCT_HEATMAP {

    input:
    path( ch_final_annots, stageAs: "raw-annotations.tsv")
    val(groupby_column)
    //Placeholder for a single product script
    path( ch_make_product_script)

    output:
    path( "product.html" ), emit: product_html
    path( "product.tsv" ), emit: product_tsv

    script:
    """
    python -m ${ch_make_product_script} \\
    --annotations ${ch_final_annots} \\
    --groupby-column ${groupby_column} \\

    """
}
