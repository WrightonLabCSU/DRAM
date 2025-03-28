process PARSE_HMM {
    
    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"

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