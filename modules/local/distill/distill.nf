process DISTILL {
    label 'process_medium'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_openpyxl_click:71fdc06a3fdfbdc1"

    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )
    path( ch_rrna_combined, stageAs: "rrna_combined.tsv" )
    path( ch_trna_combined, stageAs: "trna_combined.tsv" )

    output:
    path( "distillate.xlsx" ), emit: distillate
    path( "dram.log" ), emit: distill_log
    path( "summarized_genomes.tsv" ), emit: summarized_genomes

    script:
    """
    # export constants for script
    export FASTA_COLUMN="${params.CONSTANTS.FASTA_COLUMN}"

    distill.py -i ${ch_combined_annotations} --rrna_path ${ch_rrna_combined} --trna_path ${ch_trna_combined} --distil_topics "${params.distill_topic}" --distil_ecosystem "${params.distill_ecosystem}" --custom_distillate "${params.distill_custom}"

    """
}
