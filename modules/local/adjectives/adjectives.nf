process ADJECTIVES {
    label 'process_small'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_click_graphviz_pruned:9cc2e19fdda0e3ed"

    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )
    path( rules_tsv )

    output:
    path("adjectives.tsv"), emit: adjectives_ch

    script:
    def args = task.ext.args ?: ""
    def traits_ls = params.adjectives_list ? "--adjectives_list '${params.adjectives_list}'" : ""

    """
    # export constants for script
    export FASTA_COLUMN="${params.CONSTANTS.FASTA_COLUMN}"

    adjectives.py --annotations ${ch_combined_annotations} --output adjectives.tsv ${traits_ls} --rules_tsv '${rules_tsv}' ${args}
    """
}
