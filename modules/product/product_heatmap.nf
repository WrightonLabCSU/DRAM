/*
This is purely a placeholder process for future Product development and can be modified freely
*/
process PRODUCT_HEATMAP {

    input:
    path( ch_final_annots, stageAs: "raw-annotations.tsv")
    path( ch_distillate, stageAs: "distillate.xlsx")
    path( ch_etc_module_form, stageAs: "etc_module_database.tsv" )
    path( ch_function_heatmap_form, stageAs: "unction_heatmap_form.tsv" )
    path( ch_module_step_form, stageAs: "module_step_form.tsv" )
    //Placeholder for a single product script
    path( ch_make_product_script)
    //Place holder for a directory of product scripts
    path( ch_product_scripts)

    output:
    path( "product.html" ), emit: product_html
    path( "product.tsv" ), emit: product_tsv

    script:
    """
    #Symlink all the files in the scripts dir to the current dir
    #Before this is finalized, need to ensure we either bring in the forms with this or make usre they are not in this dir
    #    This will prevent issues with files with the same names.
    ln -s ${ch_product_scripts}/* .

    # Create a log directory if it doesn't exist
    mkdir -p logs

    # Define the log file path
    log_file="logs/product.log"

    python ${ch_make_product_script} \\
    --input-target-counts ${target_id_counts} \\
    --input-etc ${params.etc_mdoule_database} \\
    --input-module-step ${params.module_step_form} \\
    --input-function-heatmap ${params.function_heatmap_form} \\
    --output-file "product.tsv"  >> "\$log_file" 2>&1

    """
}
