process CANTHYD_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file )
    val( top_hit )
    file( ch_canthd_list )
    file( ch_canthyd_formatter )

    output:
    tuple val( sample ), path ( "${sample}_formatted_canthyd_hits.out" ), emit: canthyd_formatted_hits


    script:
    """
    python ${ch_canthyd_formatter} --hits_csv ${hits_file} --ch_canthyd_ko ${ch_canthyd_list} --output "${sample}_formatted_canthyd_hits.out"
    
    """
}

