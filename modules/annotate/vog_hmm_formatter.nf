process VOG_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file )
    val( top_hit )
    file( ch_vog_list )
    file( ch_vog_formatter )

    output:
    tuple val( sample ), path ( "${sample}_formatted_vog_hits.out" ), emit: vog_formatted_hits


    script:
    """
    python ${ch_vog_formatter} --hits_csv ${hits_file} --ch_vog_ko ${ch_vog_list} --output "${sample}_formatted_vog_hits.out"
    
    """
}

