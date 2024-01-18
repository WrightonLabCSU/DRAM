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
    collected_data = pd.DataFrame(columns=["gene_description", "module", "header", "subheader"] + samples)

    # Iterate through each TSV file
    for file in tsv_files:
        # Read the TSV file into a DataFrame
        df = pd.read_csv(file, sep='\t', skiprows=[1])

        # Print the columns of the DataFrame to check its structure
        print(f"Columns of DataFrame from {file}:\n{df.columns}")

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

        # Count occurrences of each gene_id for the current sample
        gene_id_counts = df['gene_id'].value_counts().reset_index()
        gene_id_counts.columns = ['gene_id', sample_name]

        # Merge collected_data with the counts
        collected_data = pd.merge(collected_data, gene_id_counts, on='gene_id', how='outer').fillna(0)

    # Reorder columns
    collected_data = collected_data[["gene_description", "module", "header", "subheader"] + samples]

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)


    """
}
