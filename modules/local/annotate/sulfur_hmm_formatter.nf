process SULFUR_HMM_FORMATTER {
    label 'process_tiny'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:4a55e4bf58e4a06b"
    
    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )

    output:
    tuple val( input_fasta ), path ( "${input_fasta}_formatted_sulfur_hits.csv" ), emit: sulfur_formatted_hits


    script:
    """
    sulfur_hmm_formatter.py --hits_csv ${hits_file} --gene_locs ${prodigal_locs_tsv} --output "${input_fasta}_formatted_sulfur_hits.csv"
    
    """
}
