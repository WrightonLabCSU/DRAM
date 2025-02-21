process PARSE_HMM {
    
    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"

    tag{ sample }
    
    input:
    tuple val( sample ), path( inputHMMSearch )

    output:
    tuple val( sample ), path( "${sample}_parsed_hmm.out" ), emit: parsed_hmm

    script:
    """
    parse_hmmsearch.py ${inputHMMSearch} ${sample}_parsed_hmm.out
    """
}