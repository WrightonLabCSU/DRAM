process KEGG_HMM_FORMATTER {

    tag { sample }

    input:
    file( hmm_info_path )
    tuple val( sample ), path( hits_file )
    val( top_hit )

    output:
    tuple val( sample ), path ( "${sample}_formatted_hits.out" ), emit: formatted_hits

    script:
    """
    python /home/rwoyda/Projects/DRAM2-Nextflow/DRAM2-NF/assets/kegg_hmm_formatter.py --hits_csv ${hits_file} --hmm_info_path ${hmm_info_path} --top_hit "${top_hit}" --output "${sample}_formatted_hits.out"
    """
}
