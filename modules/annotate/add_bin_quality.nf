process ADD_BIN_QUALITY {

    errorStrategy 'finish'

    input:
    file( combined_annotations )

    output:
    path("annots_bin_quality.tsv"), emit: annots_bin_quality_out, optional: true

    script:
    """
    #!/usr/bin/env python
    """
}
