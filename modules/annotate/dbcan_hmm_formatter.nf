process DBCAN_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file )
    val( top_hit )
    file( ch_dbcan_fam )
    file( ch_dbcan_subfam )
    file( ch_dbcan_formatter )

    output:
    tuple val( sample ), path ( "${sample}_formatted_dbcan_hits.out" ), emit: dbcan_formatted_hits


    script:
    """
    python ${ch_dbcan_formatter} --hits_csv ${hits_file} --fam ${ch_dbcan_fam} --subfam ${ch_dbcan_subfam} --output "${sample}_formatted_dbcan_hits.out"
    
    """
}
