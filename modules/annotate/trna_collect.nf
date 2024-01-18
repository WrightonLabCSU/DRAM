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
        unique_gene_ids = df['gene_id'].unique()
        counts = df['gene_id'].value_counts()

        for gene_id in unique_gene_ids:
            count = counts[gene_id] if gene_id in counts else 0
            collected_data.loc[:, 'gene_id'] = gene_id
            collected_data.loc[:, 'gene_description'] = df.loc[df['gene_id'] == gene_id, 'gene_description'].values[0]
            collected_data.loc[:, 'module'] = df.loc[df['gene_id'] == gene_id, 'module'].values[0]
            collected_data.loc[:, 'header'] = df.loc[df['gene_id'] == gene_id, 'header'].values[0]
            collected_data.loc[:, 'subheader'] = df.loc[df['gene_id'] == gene_id, 'subheader'].values[0]
            collected_data.loc[:, sample_name] = count

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)
    """
}
