process PRODUCT_HEATMAP {

    input:
    path( target_id_counts )

    output:
    //path( "product.html" ), emit: product_html
    path( "product.tsv" ), emit: product_tsv

    script:
    """
    # Create a log directory if it doesn't exist
    mkdir -p logs
  
    # Define the log file path
    log_file="logs/combine_annotations.log"

    python /home/rwoyda/Projects/DRAM2-Nextflow/DRAM2-NF/assets/product.py \\
    --input-target-counts ${target_id_counts} \\
    --input-etc ${params.etc_mdoule_database} \\
    --input-module-step ${params.module_step_form} \\
    --input-function-heatmap ${params.function_heatmap_form} \\
    --output-file "product.tsv"  >> "\$log_file" 2>&1

    """
}
