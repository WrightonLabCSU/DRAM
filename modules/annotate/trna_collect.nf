process TRNA_COLLECT {

    errorStrategy 'finish'

    input:
    file combined_trnas

    output:
    tuple val(sample), path("collected_trnas.tsv"), emit: trna_collected_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    def extract_samples_and_paths(combined_trnas):
        samples_and_paths = []
        for i in range(0, len(combined_trnas), 2):
            sample = combined_trnas[i].strip('[], ')
            path = combined_trnas[i + 1].strip('[], ')
            samples_and_paths.append((sample, path))
        return samples_and_paths

    # Extract sample names and paths from combined_trnas
    combined_trnas_list = extract_samples_and_paths("${combined_trnas}")

    # Create an empty DataFrame to store the collected data
    collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "module", "header", "subheader"] + [sample for sample, _ in combined_trnas_list])

    # Iterate over each sample and corresponding path
    for sample, path in combined_trnas_list:

        # Read the processed tRNAs file for the current sample
        try:
            trna_data = pd.read_csv(path, sep="\\t", skiprows=[0, 2])
        except pd.errors.EmptyDataError:
            continue  # Skip empty files

        # Add data to the collected DataFrame
        # Update the following line based on how you want to populate the values
        # collected_data[sample] = ...

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\\t", index=False)
    """
}
