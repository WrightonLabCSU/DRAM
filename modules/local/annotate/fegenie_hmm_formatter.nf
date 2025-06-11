process FEGENIE_HMM_FORMATTER {
    label 'process_tiny'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:fc59f737a5e0566a"
    
    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )

    output:
    tuple val( input_fasta ), path ( "${input_fasta}_formatted_fegenie_hits.csv" ), emit: fegenie_formatted_hits


    script:
    """
    fegenie_hmm_formatter.py --hits_csv ${hits_file} --gene_locs ${prodigal_locs_tsv} --output "${input_fasta}_formatted_fegenie_hits.csv"
    
    """
}
