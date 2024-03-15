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

    # Initialize a flag to check if all files are "NULL"
    all_files_null = True

    # Extract sample names from the file names
    samples = [os.path.basename(file).replace("_processed_trnas.tsv", "") for file in tsv_files]

    # Create an empty DataFrame to store the collected data
    collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "category", "topic_ecosystem", "subcategory"] + samples)

    # Create a dictionary to store counts for each sample
    sample_counts = {sample: [] for sample in samples}

    # Iterate through each input file
    for file in tsv_files:
        with open(file, 'r') as f:
            contents = f.read().strip()
            # Check if the file content is "NULL"
            if contents == "NULL":
                continue  # Skip this file
            all_files_null = False  # At least one file is not "NULL"
            # Read the input file into a DataFrame
            input_data = pd.read_csv(file, sep='\t')

            # Populate the gene_id column
            collected_data = pd.concat([collected_data, input_data[['gene_id']].drop_duplicates()], ignore_index=True)

            # Count occurrences for each sample
            for sample in samples:
                sample_counts[sample].extend(input_data[input_data['sample'] == sample]['gene_id'].tolist())

    # Check if all files are "NULL"
    if all_files_null:
        with open("collected_trnas.tsv", "w") as f:
            f.write("NULL")
    else:
        # Continue processing
        collected_data['type'] = ""
        collected_data['codon'] = ""
        collected_data['gene_description'] = collected_data['gene_id'].apply(lambda x: x.split('(')[0].strip() + " tRNA")
        collected_data['category'] = "tRNA"
        collected_data['topic_ecosystem'] = "tRNA"
        collected_data['subcategory'] = ""
        collected_data.drop_duplicates(subset=['gene_id'], inplace=True)

        for sample in samples:
            collected_data[sample] = collected_data['gene_id'].map(lambda x: sample_counts[sample].count(x) if x in sample_counts[sample] else 0)

        collected_data.drop(['type', 'codon'], axis=1, inplace=True)
        collected_data = collected_data[collected_data['gene_id'] != 'gene_id']
        collected_data.sort_values(by='gene_id', inplace=True)
        collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)

    """
}
