process ADJECTIVES {
    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )
    path( ch_adjectives_script )

    output:
    path "sample_adjectives.tsv" into adjectives_ch

    script:
    """
    python ${ch_adjectives_script} --annotations ${ch_combined_annotations}

    """
}