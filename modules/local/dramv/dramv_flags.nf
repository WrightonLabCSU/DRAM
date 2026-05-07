process DRAMV_FLAGS {
    label 'process_small'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_click_polars_pyarrow_pruned:00822989eabb8b47"

    input:
    path annotations
    path(fastas, stageAs: "fastas/*")
    path vog_list
    path(genomad_genes, stageAs: "genomad/*")  // 0+ files; assets/NO_FILE sentinel when absent

    output:
    path "annotations_with_flags.tsv", emit: combined_annotations_with_flags
    path "*.log", emit: log, optional: true

    script:
    def length_from_end = params.amg_length_from_end ?: 5000
    // genomad_genes is staged into ./genomad/. Pass each non-sentinel file with --genomad_genes.
    def genomad_files = (genomad_genes instanceof List ? genomad_genes : [genomad_genes])
                            .findAll { it.name != 'NO_FILE' }
    def genomad_args = genomad_files.collect { "--genomad_genes genomad/${it.name}" }.join(' ')
    """
    cat fastas/* > _catalog.fa

    dramv_flags.py \\
        -i ${annotations} \\
        -o annotations_with_flags.tsv \\
        --catalog_fasta _catalog.fa \\
        --length_from_end ${length_from_end} \\
        --vog_list ${vog_list} \\
        ${genomad_args}
    """
}
