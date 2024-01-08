process KOFAM_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file )
    val( top_hit )
    file( ch_kofam_fam )
    file( ch_kofam_subfam )
    file( ch_kofam_formatter )

    output:
    tuple val( sample ), path ( "${sample}_formatted_kofam_hits.out" ), emit: kofam_formatted_hits


    script:
    """
    python ${ch_kofam_formatter} --hits_csv ${hits_file} --fam ${ch_kofam_fam} --subfam ${ch_kofam_subfam} --output "${sample}_formatted_kofam_hits.out"
    
    """
}

