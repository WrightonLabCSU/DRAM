/*
This is purely a placeholder process for future Product development and can be modified freely
*/
process PRODUCT_HEATMAP {

    input:
    path( ch_final_annots, stageAs: "raw-annotations.tsv")
    //path( ch_etc_module_form, stageAs: "etc_module_database.tsv" )
    //path( ch_function_heatmap_form, stageAs: "function_heatmap_form.tsv" )
    //path( ch_module_step_form, stageAs: "module_step_form.tsv" )
    val(groupby_column)
    //Placeholder for a single product script
    path( ch_make_product_script)

    output:
    path( "product.html" ), emit: product_html
    path( "product.tsv" ), emit: product_tsv

    script:
    """
    #Symlink all the files in the scripts dir to the current dir
    #Before this is finalized, need to ensure we either bring in the forms with this or make usre they are not in this dir
    #    This will prevent issues with files with the same names.

    python ${ch_make_product_script} \\
    --annotations ${ch_final_annots} \\
    --groupby-column ${groupby_column} \\
    --output-dir ${params.outdir}

    """
}
//     --module_steps_form ${ch_module_step_form} \\
//     --etc_steps_form ${ch_etc_module_form} \\
//     --function_steps_form ${ch_function_heatmap_form}
