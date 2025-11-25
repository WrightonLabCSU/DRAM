process SUMMARIZE {
    label 'process_medium'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_openpyxl_click_dram-viz:bd6f4fb065d73a68"

    input:
    path( combined_annotations, stageAs: "raw-annotations.tsv" )
    path( rrna_collected, stageAs: "rrna_combined.tsv" )
    path( trna_collected, stageAs: "trna_combined.tsv" )
    path( quast_stats )
    val( distill_topic )
    val( distill_ecosystem )
    val( distill_custom )

    output:
    path( "metabolism_summary.xlsx" ), emit: distillate
    path( "*.log" ), emit: log
    path( "summarized_genomes.tsv" ), emit: summarized_genomes
    path( "genome_stats.tsv" ), emit: genome_stats

    script:
    def args = task.ext.args ?: ""
    def rrna = rrna_collected ? "--rrna_path ${rrna_collected}" : ""
    def trna = trna_collected ? "--trna_path ${trna_collected}" : ""
    def quast = quast_stats ? "--trna_path ${quast_stats}" : ""

    """
    distill.py \
        -i ${combined_annotations} \
        ${rrna} \
        ${trna} \
        ${quast} \
        --distil_topics "${distill_topic}" \
        --distil_ecosystem "${distill_ecosystem}" \
        --custom_distillate "${distill_custom}" \
        ${args}
    """
}
