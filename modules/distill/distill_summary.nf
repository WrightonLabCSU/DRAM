process DISTILL_SUMMARY {

    input:
    file( combined_annotations )
    path( distill_sheets_file )
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
        --distill_sheets "${distill_sheets_file}" \
        --output "genome_summary.tsv" >> "\$log_file" 2>&1
    """
}






/*
    python "${ch_distill_summary_script}" \
        --combined_annotations "${combined_annotations}" \
        --distill_sheets "${distill_sheets}" \
        --output "genome_summary.tsv" >> "\$log_file" 2>&1


*/

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
  

    python "${ch_distill_summary_script}" \
        --combined_annotations "${combined_annotations}" \
        --genome_summary_form "${params.genome_summary_form}" \
        --output "genome_summary.tsv" >> "\$log_file" 2>&1

    """
}


*/