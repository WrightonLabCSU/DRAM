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

        # Update the count for each gene_id and sample
        for index, row in df.iterrows():
            gene_id = row['gene_id']

            # Check if the gene_id is already present in collected_data
            if gene_id in collected_data['gene_id'].values:
                # Update the count for the corresponding sample
                collected_data.loc[collected_data['gene_id'] == gene_id, sample_name] += 1
            else:
                # Gene_id not present, add a new row to collected_data
                new_row = pd.DataFrame({
                    'gene_id': [gene_id],
                    'gene_description': [row['gene_description']],
                    'module': [row['module']],
                    'header': [row['header']],
                    'subheader': [row['subheader']],
                    sample_name: [1]
                })
                collected_data = pd.concat([collected_data, new_row], ignore_index=True)

    # Fill NaN values with 0 in the sample columns
    collected_data[samples] = collected_data[samples].fillna(0).astype(int)

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)

    """
}
