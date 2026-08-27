process PRODUCT_HEATMAP {
    label 'process_small'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'oras://community.wave.seqera.io/library/python_dram-viz:c86af82699067899' :
        'community.wave.seqera.io/library/python_dram-viz:fad89247bbf19865' }"

    input:
    path(ch_final_annots, stageAs: "raw-annotations.tsv")
    val(fasta_column)
    path(rules_tsv)
    path(mapping_file)
    val(rules_system)  // comma seperated list

    output:
    path( "*.html" ), emit: product_html

    script:
    def args = task.ext.args ?: ''
    def viz_rules_tsv = rules_tsv ? "--rules_tsv $rules_tsv" : ''
    def viz_mapping_file = mapping_file ? "--mapping $mapping_file" : ''

    def rules_system_values = rules_system?.trim() ?
        rules_system.split(',').collect { it.trim() } :
        ['']

    def dram_viz_commands = rules_system_values.collect { rules_system_value ->
        def viz_rules_system = rules_system_value ? "--rules_system $rules_system_value" : ''
        def output_prefix = rules_system_value ? "product_${rules_system_value}_" : 'product_'
        """
        dram_viz \\
            --annotations ${ch_final_annots} \\
            --fasta_column ${fasta_column} \\
            $viz_rules_tsv \\
            $viz_mapping_file \\
            $viz_rules_system \\
            $args

        for product_file in product_*.html; do
            mv "\$product_file" "product_outputs/${output_prefix}\${product_file#product_}"
        done
        """
    }.join("\n")

    """
    mkdir -p product_outputs
    $dram_viz_commands
    mv product_outputs/*.html .
    """
}
