process PARSE_HMM {
    label 'process_tiny'
    
    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:4a55e4bf58e4a06b"

    tag{ input_fasta }
    
    input:
    tuple val( input_fasta ), path( inputHMMSearch )

    output:
    tuple val( input_fasta ), path( "${input_fasta}_parsed_hmm.out" ), emit: parsed_hmm

    script:
    """
    parse_hmmsearch.py ${inputHMMSearch} ${input_fasta}_parsed_hmm.out
    """
}