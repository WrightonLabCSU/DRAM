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
    tsv_files = [f for f in os.listdir('.') if f.endswith('_processed_trnas.tsv')]

    # Create an empty DataFrame to store the collected data
    collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "module", "header", "subheader"])

    # Iterate through each TSV file
    for file in tsv_files:
        # Read the TSV file into a DataFrame
        df = pd.read_csv(file, sep='\t')

        # Extract sample name from the file name
        sample_name = os.path.basename(file).replace("_processed_trnas.tsv", "")

        # Update the count for each gene_id and sample
        for index, row in df.iterrows():
            # Construct the gene_id value to match
            gene_id_to_match = f"{row['type']} ({row['codon']})"
            
            # Check if the gene_id is already present in collected_data
            if gene_id_to_match in collected_data['gene_id'].values:
                # Update the count for the corresponding sample
                collected_data.at[collected_data['gene_id'] == gene_id_to_match, sample_name] += 1
            else:
                # Gene_id not present, add a new row to collected_data
                new_row = pd.DataFrame({
                    'gene_id': [gene_id_to_match],
                    'gene_description': [f"{row['type']} tRNA with {row['codon']} Codon"],
                    'module': [f"{row['type']} tRNA"],
                    'header': ['tRNA'],
                    'subheader': [''],
                    sample_name: [1]
                })
                collected_data = pd.concat([collected_data, new_row], ignore_index=True)

    # Fill NaN values with 0 in the sample columns
    collected_data[sample_name] = collected_data[sample_name].fillna(0).astype(int)

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)
    """
}
