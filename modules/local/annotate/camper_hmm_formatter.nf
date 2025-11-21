process CAMPER_HMM_FORMATTER {
    label 'process_tiny'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:d2c88b719ab1322c"

    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    file( ch_camper_list )

    output:
    tuple val( input_fasta ), path ( "${input_fasta}___formatted_camper_hits.csv" ), emit: camper_formatted_hits


    script:
    """
    camper_hmm_formatter.py --hits_csv ${hits_file} --ch_camper_list ${ch_camper_list} --gene_locs ${prodigal_locs_tsv} --output "${input_fasta}___formatted_camper_hits.csv"
    
    """
}
