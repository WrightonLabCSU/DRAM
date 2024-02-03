process DISTILL_SUMMARY {

    input:
    file( combined_annotations )
    path( ch_combined_distill_sheets )
    file( target_id_counts )
    file( ch_distill_summary_script )
    

    output:
    path("genome_summary.tsv"), emit: ch_genome_sum_simple

    script:
    """
    # Create a log directory if it doesn't exist
    mkdir -p logs

    # Define the log file path
    log_file="logs/distill.log"

    python "${ch_distill_summary_script}" \
        --combined_annotations "${combined_annotations}" \
        --target_id_counts "${target_id_counts}" \
        --output "genome_summary.tsv" >> "\$log_file" 2>&1

    """
}

