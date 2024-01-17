process TRNA_COLLECT {

    errorStrategy 'finish'

    input:
    val combined_trnas

    output:
    tuple val(sample), path("collected_trnas.tsv"), emit: trna_collected_out, optional: true

    script:
    """
    #!/usr/bin/env python

    # Replace single quotes with double quotes
    combined_trnas = ${combined_trnas.replace("'", '"')}

    import pandas as pd

    # Create an empty DataFrame to store collected data
    collected_df = pd.DataFrame()

    # Iterate through combined_trnas to extract sample names and paths
    for i in range(0, len(combined_trnas), 2):
        sample_name = combined_trnas[i]
        file_path = combined_trnas[i + 1]

        # Read the input file into a DataFrame
        trna_frame = pd.read_csv(file_path, sep="\\t", skiprows=[0, 2])

        # Extract relevant columns from the input DataFrame
        trna_frame = trna_frame[["sample", "query_id", "tRNA #", "begin", "end", "type", "codon", "score"]]

        # Rename columns based on sample_name
        columns_mapping = {col: f"{sample_name}_{col}" for col in trna_frame.columns[1:]}
        trna_frame.rename(columns=columns_mapping, inplace=True)

        # Merge the extracted data into the collected DataFrame
        collected_df = pd.merge(collected_df, trna_frame, how='outer', on="sample")

    # Write the collected DataFrame to the output file
    collected_df.to_csv("collected_trnas.tsv", sep="\\t", index=False)
    """
}
