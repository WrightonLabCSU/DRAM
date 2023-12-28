process COMBINE_ANNOTATIONS {

    input:
    val all_annotations

    output:
    path "combined_annotations.tsv", emit: combined_annotations_out

    script:
    """
    # Create a log directory if it doesn't exist
    mkdir logs

    # Define the log file path
    touch logs/combine_annotations.log
    log_file="logs/combine_annotations.log"

    python /home/rwoyda/Projects/DRAM2-Nextflow/DRAM2-NF/assets/combine_annotations.py --annotations ${all_annotations} --output "combined_annotations.tsv" >> \$log_file 2>&1

    """
}


/*

    # If the user provided additional annotations using --add_annotations,
    # they want to add them to the current run - combine them.
    if [[ "${params.add_annotations}" != "empty" ]]; then
        python /home/rwoyda/Projects/DRAM2-Nextflow/DRAM2-NF/assets/combine_annotations.py --annotations ${all_annotations} ${params.add_annotations} --output "combined_annotations.tsv" >> \$log_file 2>&1
    else
        python /home/rwoyda/Projects/DRAM2-Nextflow/DRAM2-NF/assets/combine_annotations.py --annotations ${all_annotations} --output "combined_annotations.tsv" >> \$log_file 2>&1
    fi
*/