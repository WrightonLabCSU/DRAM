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

    # List all tsv files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Extract sample names from the file names
    samples = [os.path.basename(file).replace("_processed_trnas.tsv", "") for file in tsv_files]

    # Create an empty DataFrame to store the collected data
    collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "module", "header", "subheader"] + samples)

    # Iterate through each TSV file
    for file in tsv_files:
        # Read the TSV file into a DataFrame
        df = pd.read_csv(file, sep='\t')

        # Extract sample name from the file name
        sample_name = os.path.basename(file).replace("_processed_trnas.tsv", "")

        # Create a DataFrame to count occurrences of each unique gene_id
        gene_counts = df.groupby(['type', 'codon']).size().reset_index(name='count')

        # Merge gene counts into the main collected_data DataFrame
        collected_data = pd.merge(collected_data, gene_counts, how='left', left_on=['gene_id'], right_on=['type', 'codon'])

        # Rename the count column to the sample_name
        collected_data.rename(columns={'count': sample_name}, inplace=True)

        # Drop the unnecessary columns from the merged DataFrame
        collected_data.drop(['type', 'codon'], axis=1, inplace=True)

    # Fill NaN values with 0 in the sample columns
    collected_data[samples] = collected_data[samples].fillna(0).astype(int)

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)

    """
}
