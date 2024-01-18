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

        # Count occurrences of each unique gene_id value in the current file
        gene_id_counts = df['gene_id'].value_counts()

        # Iterate through unique gene_ids in the current file and update collected_data
        for unique_gene_id, count in gene_id_counts.items():
            # Check if the gene_id already exists in collected_data
            mask = collected_data['gene_id'] == unique_gene_id

            if not mask.any():
                # If gene_id is not in collected_data, add a new row
                new_row = pd.DataFrame([[unique_gene_id, df.loc[df['gene_id'] == unique_gene_id, 'gene_description'].values[0], df.loc[df['gene_id'] == unique_gene_id, 'module'].values[0], df.loc[df['gene_id'] == unique_gene_id, 'header'].values[0], df.loc[df['gene_id'] == unique_gene_id, 'subheader'].values[0]] + [0] * len(samples)], columns=collected_data.columns)
                collected_data = pd.concat([collected_data, new_row], ignore_index=True)

            # Update the count for the corresponding sample-named column using the index
            collected_data.at[mask[mask].index[0], sample_name] += count

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)
    """
}
