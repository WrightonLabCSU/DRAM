process GENERIC_HMM_FORMATTER {
    
    input:
    path( hmm_info_path )
    path( hits_file )
    bool( top_hit )
    
    output:
    path ( "formatted_hits.out" ), emit: formatted_hits
    
    script:
    """
    python /home/rwoyda/Projects/DRAM2-Nextflow/DRAM2-NF/assets/generic_hmm_formatter.py \
        --hits_csv ${hits_csv} \
        --hmm_info_path ${hmm_info_path} \
        --top_hit ${top_hit} \
        --output "formatted_hits.out"
    """
}