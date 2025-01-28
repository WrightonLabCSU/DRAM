process PARSE_HMM {

    tag{ sample }
    
    input:
    tuple val( sample ), path( inputHMMSearch )
    file( ch_parse_hmmsearch )

    output:
    tuple val( sample ), path( "${sample}_parsed_hmm.out" ), emit: parsed_hmm

    script:
    """
    python ./${ch_parse_hmmsearch} ${inputHMMSearch} ${sample}_parsed_hmm.out
    """
}