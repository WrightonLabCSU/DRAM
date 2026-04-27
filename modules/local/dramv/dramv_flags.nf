process DRAMV_FLAGS {
    label 'process_small'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_click_polars_pyarrow_pruned:00822989eabb8b47"

    input:
    path annotations
    path(fastas, stageAs: "fastas/*")

    output:
    path "annotations_with_flags.tsv", emit: combined_annotations_with_flags
    path "*.log", emit: log, optional: true

    script:
    def length_from_end = params.amg_length_from_end ?: 5000
    """
    cat fastas/* > _catalog.fa

    dramv_flags.py \\
        -i ${annotations} \\
        -o annotations_with_flags.tsv \\
        --catalog_fasta _catalog.fa \\
        --length_from_end ${length_from_end}
    """
}
