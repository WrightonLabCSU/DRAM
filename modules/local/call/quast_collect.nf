process QUAST_COLLECT {

    errorStrategy 'finish'

    input:
    file quast_tsv_files

    output:
    path("collected_quast.tsv"), emit: quast_collected_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import os
    import pandas as pd
    import re

    # List all tsv files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Extract sample names from the file names
    samples = [os.path.basename(file).replace("_report.tsv", "") for file in tsv_files]


    """

}