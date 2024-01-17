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

    # Extract sample names and paths from combined_trnas
    combined_trnas_list = ${combined_trnas.toList()}
    samples_and_paths = [combined_trnas_list[i:i+2] for i in range(0, len(combined_trnas_list), 2)]

    # Create an empty DataFrame to store the collected data
    collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "module", "header", "subheader"] + [sample[0] for sample in samples_and_paths])

    # Iterate over each sample and corresponding path
    for sample, path in samples_and_paths:
        # Read the processed tRNAs file for the current sample
        try:
            trna_data = pd.read_csv(path, sep="\t", skiprows=[0, 2])
        except pd.errors.EmptyDataError:
            continue  # Skip empty files

        # Add data to the collected DataFrame
        # Update the following line based on how you want to populate the values
        # collected_data[sample] = ...

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)
    """
}
