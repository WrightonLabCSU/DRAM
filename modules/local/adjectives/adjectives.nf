process ADJECTIVES {
    label 'process_low'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_click_graphviz_pruned:34469eb76a8384ac"

    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )

    output:
    path("adjectives.tsv"), emit: adjectives_ch

    script:
    """
    # export constants for script
    export FASTA_COLUMN="${params.CONSTANTS.FASTA_COLUMN}"

    adjectives.py --annotations ${ch_combined_annotations} --output adjectives.tsv --adjectives_list '${params.adjectives_list}'

    """
}
