process ADJECTIVES {
    label 'process_low'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_click_graphviz_pruned:34469eb76a8384ac"

    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )
    path( ch_adjectives_script )

    output:
    path "input_fasta_adjectives.tsv" into adjectives_ch

    script:
    """
    rule_adjectives --annotations ${ch_combined_annotations}

    """
}
