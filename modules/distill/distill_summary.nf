process DISTILL_SUMMARY {
    input:
    file( combined_annotations )
    file( target_id_counts )
    file( ch_distill_summary_script )

    output:
    path("genome_summary.tsv"), emit: metab_summ_simple

    script:
    """
    # Create a log directory if it doesn't exist
    mkdir -p logs

    # Define the log file path
    log_file="logs/distill.log"

    python "${ch_distill_summary_script}" \
        --combined_annotations "${combined_annotations}" \
        --genome_summary_form "${params.genome_summary_form}" \
        --output "genome_summary.tsv" >> "\$log_file" 2>&1

    """
}

/*
process DISTILL_SUMMARY {

    input:
    file( combined_annotations )
    file( target_id_counts )
    file( ch_distill_summary_script )

    output:
    path("genome_summary.tsv"), emit: metab_summ_simple

    script:
    """
    # Create a log directory if it doesn't exist
    mkdir -p logs
  
    # Define the log file path
    log_file="logs/distill.log"
  
    if [[ "${params.add_module1}" != "empty" || "${params.add_module2}" != "empty" || "${params.add_module3}" != "empty" || "${params.add_module4}" != "empty" || "${params.add_module5}" != "empty" ]]; then
        python "${ch_distill_summary_script}" \
            --combined_annotations "${combined_annotations}" \
            --genome_summary_form "${params.genome_summary_form}" \
            --target_id_counts "${target_id_counts}" \
            --output "genome_summary.tsv" \
            --add_module1 "${params.add_module1}" \
            --add_module2 "${params.add_module2}" \
            --add_module3 "${params.add_module3}" \
            --add_module4 "${params.add_module4}" \
            --add_module5 "${params.add_module5}" >> "\$log_file" 2>&1
    else
        python "${ch_distill_summary_script}" \
            --combined_annotations "${combined_annotations}" \
            --genome_summary_form "${params.genome_summary_form}" \
            --target_id_counts "${target_id_counts}" \
            --output "genome_summary.tsv" >> "\$log_file" 2>&1
    fi
    """
}
*/