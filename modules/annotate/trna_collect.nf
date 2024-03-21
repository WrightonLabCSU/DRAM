#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TRNA_COLLECT {

    errorStrategy 'finish'

    input:
    file combined_trnas

    output:
    path("collected_trnas.tsv"), emit: trna_collected_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import os
    import pandas as pd
    import re

    # List all tsv files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Check if all file contents are "NULL"
    all_files_null = True
    for file in tsv_files:
        with open(file, 'r') as f:
            contents = f.read().strip()
            if contents != "NULL":
                all_files_null = False
                break

    # If all file contents are "NULL", set the output file's content to "NULL"
    if all_files_null:
        with open('collected_trnas.tsv', 'w') as f:
            f.write("NULL")
    else:
        # Process files as usual if not all contents are "NULL"

        # Extract sample names from the file names
        samples = [os.path.basename(file).replace("_processed_trnas.tsv", "") for file in tsv_files]

        # Create an empty DataFrame to store the collected data
        collected_data = pd.DataFrame()

        # Iterate through each input file
        for file in tsv_files:
            # Read the input file into a DataFrame, if not "NULL"
            with open(file, 'r') as f:
                contents = f.read().strip()
                if contents == "NULL":
                    continue
                input_data = pd.read_csv(file, sep='\t', header=None, names=["sample", "query_id", "tRNA #", "begin", "end", "type", "codon", "score", "gene_id"])
                # Proceed with processing as before

        # Note: The rest of your processing code remains the same as in your original script,
        # including collecting data into 'collected_data' and writing it to 'collected_trnas.tsv'.
        # Make sure to adjust the logic accordingly if there are parts that depend on the initial structure of 'collected_data'.
    """
}
