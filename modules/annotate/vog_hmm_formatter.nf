process VOG_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file )
    val( top_hit )
    file( ch_vog_fam )
    file( ch_vog_subfam )
    file( ch_vog_formatter )

    output:
    tuple val( sample ), path ( "${sample}_formatted_vog_hits.out" ), emit: vog_formatted_hits


    script:
    """
    python ${ch_vog_formatter} --hits_csv ${hits_file} --fam ${ch_vog_fam} --subfam ${ch_vog_subfam} --output "${sample}_formatted_vog_hits.out"
    
    """
}

