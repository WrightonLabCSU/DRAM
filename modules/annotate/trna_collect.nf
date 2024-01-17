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

    # Convert the combined_trnas string to a list
    combined_trnas_list = ${combined_trnas}

    # Initialize an empty DataFrame to store the collected data
    collected_data = pd.DataFrame()

    # Iterate through the combined_trnas_list to extract sample names
    for i in range(0, len(combined_trnas_list), 2):
        sample_name = combined_trnas_list[i][1:-1]  # Remove single quotes from sample name
        collected_data[sample_name] = ""

    # Write the DataFrame to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\\t", index=False)
    """
}
