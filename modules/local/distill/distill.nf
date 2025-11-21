process SUMMARIZE {
    label 'process_medium'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_openpyxl_click_dram-viz:bd6f4fb065d73a68"

    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )
    path( ch_rrna_collected, stageAs: "rrna_combined.tsv" )
    path( ch_trna_collected, stageAs: "trna_combined.tsv" )
    path( ch_quast_stats )
    val( distill_topic )
    val( distill_ecosystem )
    val( distill_custom )

    output:
    path( "metabolism_summary.xlsx" ), emit: distillate
    path( "*.log" ), emit: log
    path( "summarized_genomes.tsv" ), emit: summarized_genomes
    path( "genome_stats.tsv" ), emit: genome_stats

    script:
    """
    # export constants for script
    export FASTA_COLUMN="${params.CONSTANTS.FASTA_COLUMN}"

    distill.py -i ${ch_combined_annotations} --rrna_path '${ch_rrna_collected}' --trna_path '${ch_trna_collected}' --distil_topics "${distill_topic}" --distil_ecosystem "${distill_ecosystem}" --custom_distillate "${distill_custom}" --quast_path '${ch_quast_stats}'

    """
}
