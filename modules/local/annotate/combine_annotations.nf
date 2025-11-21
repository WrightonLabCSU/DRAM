process COMBINE_ANNOTATIONS {
    label 'process_small'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:d2c88b719ab1322c"

    input:
    path(fastas, stageAs: "annotations/*" )
    path(genes, stageAs: "genes/*" )

    output:
    path "raw-annotations.tsv", emit: combined_annotations_out
    path( "*.log" ), emit: log

    script:
    """
    export FASTA_COLUMN="${params.CONSTANTS.FASTA_COLUMN}"

    combine_annotations.py \\
        --annotations_dir annotations \\
        --genes_dir genes \\
        --threads "${params.threads}" \\
        --output raw-annotations.tsv

    """
}
