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
    path( ch_make_product_script)

    output:
    path( "product.html" ), emit: product_html
    path( "product.tsv" ), emit: product_tsv

    script:
    """
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
