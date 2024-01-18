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

    # Create a dictionary to store counts for each sample
    sample_counts = {sample: [] for sample in samples}

    # Iterate through each input file
    for file in tsv_files:
        # Read the input file into a DataFrame
        input_data = pd.read_csv(file, sep='\t', header=None, names=["sample", "query_id", "tRNA #", "begin", "end", "type", "codon", "score", "gene_id"])

        # Populate the gene_id column
        collected_data = pd.concat([collected_data, input_data[['gene_id']].drop_duplicates()], ignore_index=True)

        # Count occurrences for each sample
        for sample in samples:
            sample_counts[sample].extend(input_data[input_data['sample'] == sample]['gene_id'])

    # Add 'type' and 'codon' columns to collected_data
    collected_data['type'] = ""
    collected_data['codon'] = ""

    # Populate other columns based on the given rules
    collected_data['gene_description'] = collected_data['gene_id'].replace(r'\([^)]*\)', '', regex=True) + " tRNA with " + collected_data['codon'] + " Codon"
    collected_data['module'] = collected_data['gene_id'].replace(r'\([^)]*\)', '', regex=True) + " tRNA"
    collected_data['header'] = "tRNA"
    collected_data['subheader'] = ""

    # Deduplicate the rows based on gene_id
    collected_data.drop_duplicates(subset=['gene_id'], inplace=True)

    # Count occurrences for each sample and fill the additional columns dynamically
    for sample in samples:
        collected_data[sample] = collected_data['gene_id'].map(lambda x: sample_counts[sample].count(x) if x in sample_counts[sample] else 0)

    # Remove 'type' and 'codon' columns
    collected_data.drop(['type', 'codon'], axis=1, inplace=True)

    # Remove the first row that includes extra column names
    collected_data = collected_data[collected_data['gene_id'] != 'gene_id']

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)
    """
}
