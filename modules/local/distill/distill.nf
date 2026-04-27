process SUMMARIZE {
    label 'process_medium'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_click_polars_pyarrow_pruned:00822989eabb8b47"

    input:
    path combined_annotations
    path rrna_collected
    path trna_collected
    path quast_stats
    val distill_topic
    val distill_ecosystem
    val distill_custom

    output:
    path "metabolism_summary.xlsx", emit: distillate
    path "*.log", emit: log
    path "summarized_genomes.tsv", emit: summarized_genomes
    path "genome_stats.tsv", emit: genome_stats

    script:
    def args = task.ext.args ?: ""
    def rrna = rrna_collected ? "--rrna_path ${rrna_collected}" : ""
    def trna = trna_collected ? "--trna_path ${trna_collected}" : ""
    def quast = quast_stats ? "--quast_path ${quast_stats}" : ""
    def groupby = params.use_dramv ? "scaffold" : params.groupby_column
    def amg_only = params.use_dramv ? "--amg_only" : ""

    """
    distill.py \
        -i ${combined_annotations} \
        ${rrna} \
        ${trna} \
        ${quast} \
        --groupby_column ${groupby} \
        --distil_topics "${distill_topic}" \
        --distil_ecosystem "${distill_ecosystem}" \
        --custom_distillate "${distill_custom}" \
        ${amg_only} \
        ${args}
    """
}
