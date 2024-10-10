process CANTHYD_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    val( top_hit )
    file( ch_canthyd_list )
    file( ch_canthyd_formatter )

    output:
    tuple val( sample ), path ( "${sample}_formatted_canthyd_hits.csv" ), emit: canthyd_formatted_hits


    script:
    """
    python ${ch_canthyd_formatter} --hits_csv ${hits_file} --ch_canthyd_ko ${ch_canthyd_list} --gene_locs ${prodigal_locs_tsv} --output "${sample}_formatted_canthyd_hits.csv"
    
    """
}

