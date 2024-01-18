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
        df['gene_id'] = df.apply(lambda row: f"{row['type']} ({row['codon']})" if 'pseudo' not in row['type'].lower() else f"{row['type']} (pseudo) ({row['codon']})", axis=1)

        # Construct other columns as before
        df['gene_description'] = df['type'] + ' tRNA with ' + df['codon'] + ' Codon'
        df['module'] = df['type'] + ' tRNA'
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
        # Leave sample-named columns empty for now
        collected_data.loc[:, sample_name] = ''

    # Extract unique variations of gene_id
    unique_gene_ids = collected_data['gene_id'].unique()

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)


    """
}
