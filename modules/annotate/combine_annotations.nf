process COMBINE_ANNOTATIONS {

    input:
    val all_annotations
    file(ch_combine_annot_script)

    output:
    path "raw-annotations.tsv", emit: combined_annotations_out

    script:
    """
    # Create a log directory if it doesn't exist
    mkdir logs

    # Define the log file path
    touch logs/combine_annotations.log
    log_file="logs/combine_annotations.log"

    python ${ch_combine_annot_script} --annotations ${all_annotations} --threads ${params.threads} --output "raw-annotations.tsv" >> \$log_file 2>&1

    """
}
