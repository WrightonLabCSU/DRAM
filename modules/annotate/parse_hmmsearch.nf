process PARSE_HMM {

    tag{ sample }
    
    input:
    tuple val( sample ), path( inputHMMSearch )

    output:
    tuple val( sample ), path( "${sample}_parsed_hmm.out" ), emit: parsed_hmm

    script:
    """
    bash -c "python \$(readlink -f ./assets/parse_hmmsearch.py) small_sample-fasta_hmmsearch.out ${sample}_parsed_hmm.out"

    """
}