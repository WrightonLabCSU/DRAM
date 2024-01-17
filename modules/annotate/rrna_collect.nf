process RRNA_COLLECT {

    errorStrategy 'finish'

    input:
    val combined_rrnas

    output:
    tuple val(sample), path("${sample}_collected_rrnas.tsv"), emit: rrna_collected_out, optional: true

    script:
    """
    #!/usr/bin/env python


    """
}
