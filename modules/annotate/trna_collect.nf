process TRNA_COLLECT {

    errorStrategy 'finish'

    input:
    val combined_trnas

    output:
    tuple val(sample), path("collected_trnas.tsv"), emit: trna_collected_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    # Create a dictionary to store data for each sample
    sample_data = {}

    # Loop through combined_trnas to extract sample names and paths
    for i in range(0, len(combined_trnas), 2):
        sample_name = combined_trnas[i].strip("'")
        sample_path = combined_trnas[i + 1]

        # Create an empty DataFrame for each sample
        sample_data[sample_name] = pd.read_csv(sample_path, sep="\\t", nrows=0)

    # Concatenate DataFrames along columns to form the collected_trnas DataFrame
    collected_trnas = pd.concat(sample_data.values(), axis=1)

    # Write the collected_trnas DataFrame to the output file
    collected_trnas.to_csv("collected_trnas.tsv", sep="\\t", index=False)
    """
}
