process PARSE_HMM {

    tag{ sample }
    
    input:
    tuple val( sample ), path( inputHMMSearch )

    output:
    tuple val( sample ), path( "${sample}_parsed_hmm.out" ), emit: parsed_hmm

    script:
    """
    python ./assets/parse_hmmsearch.py ${inputHMMSearch} ${sample}_parsed_hmm.out
    """
}