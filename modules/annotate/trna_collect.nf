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

    # Extract sample names and paths from the combined_trnas variable
    combined_trnas_list = ${combined_trnas}
    samples = [combined_trnas_list[i].strip("'") for i in range(0, len(combined_trnas_list), 2)]
    paths = [combined_trnas_list[i + 1] for i in range(0, len(combined_trnas_list), 2)]

    # Create an empty DataFrame with the desired columns
    collected_trnas = pd.DataFrame(columns=["gene_id", "gene_description", "module", "header", "subheader"] + samples)

    # Save the DataFrame to the output file
    collected_trnas.to_csv("collected_trnas.tsv", sep="\\t", index=False)
    """
}
