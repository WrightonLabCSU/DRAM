process KOFAM_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    val( top_hit )
    file( ch_kofam_list )
    file( ch_kofam_formatter )

    output:
    tuple val( sample ), path ( "${sample}_formatted_kofam_hits.csv" ), emit: kofam_formatted_hits


    script:
    """
    python ${ch_kofam_formatter} --hits_csv ${hits_file} --ch_kofam_ko ${ch_kofam_list} --gene_locs "gene_locs.tsv" --output "${sample}_formatted_kofam_hits.csv"
    
    """
}

