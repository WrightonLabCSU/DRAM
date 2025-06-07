// This isn't used in the current version, but we need to fix the absolute paths in the script block
process GENERIC_HMM_FORMATTER {
    label 'process_tiny'
    
    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:4a55e4bf58e4a06b"

    input:
    path( hmm_info_path )
    path( hits_file )
    bool( top_hit )
    
    output:
    path ( "formatted_hits.out" ), emit: formatted_hits
    
    script:
    """
    generic_hmm_formatter.py \
        --hits_csv ${hits_csv} \
        --hmm_info_path ${hmm_info_path} \
        --top_hit ${top_hit} \
        --output "formatted_hits.out"
    """
}