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

    # Create an empty DataFrame to store the collected data
    collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "module", "header", "subheader"] + [os.path.basename(file).replace("_processed_trnas.tsv", "") for file in tsv_files])

    # Global counter for gene_id occurrences
    gene_id_counter = {}

    # Iterate through each TSV file
    for file in tsv_files:
        # Read the TSV file into a DataFrame
        df = pd.read_csv(file, sep='\t', skiprows=[1])

        # Construct the gene_id column
        df['gene_id'] = df['type'] + ' (' + df['codon'] + ')'

        # Construct the gene_description column
        df['gene_description'] = df['type'] + ' tRNA with ' + df['codon'] + ' Codon'

        # Construct the module column
        df['module'] = df['type'] + ' tRNA'

        # Add constant values to header and subheader columns
        df['header'] = 'tRNA'
        df['subheader'] = ''

        # Extract sample name from the file name
        sample_name = os.path.basename(file).replace("_processed_trnas.tsv", "")

        # Update the corresponding columns in the collected_data DataFrame
        collected_data.loc[:, 'gene_id'] = df['gene_id']
        collected_data.loc[:, 'gene_description'] = df['gene_description']
        collected_data.loc[:, 'module'] = df['module']
        collected_data.loc[:, 'header'] = df['header']
        collected_data.loc[:, 'subheader'] = df['subheader']

        # Update the global gene_id_counter
        for unique_gene_id in df['gene_id'].unique():
            gene_id_counter[unique_gene_id] = gene_id_counter.get(unique_gene_id, 0) + 1

    # Update the sample-named columns based on the global gene_id_counter
    for unique_gene_id, count in gene_id_counter.items():
        collected_data.loc[collected_data['gene_id'] == unique_gene_id, tsv_files] = count

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)

    """
}
