process KEGG_HMM_FORMATTER {

    tag { sample }

    input:
    file( hmm_info_path )
    tuple val( sample ), path( hits_file )
    val( top_hit )

    output:
    tuple val( sample ), path ( "${sample}_formatted_kegg_hits.csv" ), emit: formatted_hits

    script:
    """
    kegg_hmm_formatter.py --hits_csv ${hits_file} --hmm_info_path ${hmm_info_path} --top_hit "${top_hit}" --output "${sample}_formatted_kegg_hits.csv"
    """
}
