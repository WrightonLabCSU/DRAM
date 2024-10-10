process CAMPER_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    val( top_hit )
    file( ch_camper_list )
    file( ch_camper_formatter )

    output:
    tuple val( sample ), path ( "${sample}_formatted_camper_hits.csv" ), emit: camper_formatted_hits


    script:
    """
    python ${ch_camper_formatter} --hits_csv ${hits_file} --ch_camper_list ${ch_camper_list} --gene_locs ${prodigal_locs_tsv} --output "${sample}_formatted_camper_hits.csv"
    
    """
}

