process CAMPER_HMM_FORMATTER {
    label 'process_tiny'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:7b4f27307e83be0e"

    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    file( ch_camper_list )

    output:
    tuple val( input_fasta ), path ( "${input_fasta}_formatted_camper_hits.csv" ), emit: camper_formatted_hits


    script:
    """
    camper_hmm_formatter.py --hits_csv ${hits_file} --ch_camper_list ${ch_camper_list} --gene_locs ${prodigal_locs_tsv} --output "${input_fasta}_formatted_camper_hits.csv"
    
    """
}
