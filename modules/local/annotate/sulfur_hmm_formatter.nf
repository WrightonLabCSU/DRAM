process SULFUR_HMM_FORMATTER {

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2:89f055454dac3575"
    
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
