process KOFAM_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file )
    file( ch_sulfur_formatter )

    output:
    tuple val( sample ), path ( "${sample}_formatted_sulfur_hits.csv" ), emit: sulfur_formatted_hits


    script:
    """
    python ${ch_sulfur_formatter} --hits_csv ${hits_file} --output "${sample}_formatted_sulfur_hits.csv"
    
    """
}

