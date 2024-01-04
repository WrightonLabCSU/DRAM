process PARSE_HMM {

    tag{ sample }
    
    input:
    tuple val( sample ), path( inputHMMSearch )

    output:
    tuple val( sample ), path( "${sample}_parsed_hmm.out" ), emit: parsed_hmm

    script:
    """
    python ${params.parse_hmmsearch_script} ${inputHMMSearch} ${sample}_parsed_hmm.out
    """
}