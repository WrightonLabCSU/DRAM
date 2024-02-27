process FEGENIE_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file )
    file( ch_fegenie_formatter )

    output:
    tuple val( sample ), path ( "${sample}_formatted_fegenie_hits.csv" ), emit: fegenie_formatted_hits


    script:
    """
    python ${ch_fegenie_formatter} --hits_csv ${hits_file} --output "${sample}_formatted_fegenie_hits.csv"
    
    """
}

