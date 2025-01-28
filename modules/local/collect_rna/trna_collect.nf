#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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

    # Function to check if file content is "NULL"
    def file_content_is_null(filename):
        with open(filename, 'r') as file:
            content = file.read().strip()
        return content == "NULL"

    # List all tsv files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Check if all files contain "NULL"
    all_files_null = all(file_content_is_null(file) for file in tsv_files)

    if all_files_null:
        # If all files contain "NULL", write "NULL" to the output and exit
        with open("collected_trnas.tsv", "w") as output_file:
            output_file.write("NULL")
    else:
        # Filter out files with "NULL" content
        tsv_files = [file for file in tsv_files if not file_content_is_null(file)]

        # Extract sample names from the file names
        samples = [os.path.basename(file).replace("_processed_trnas.tsv", "") for file in tsv_files]

        # Create an empty DataFrame to store the collected data
        collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "category", "topic_ecosystem", "subcategory"] + samples)

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

        # Deduplicate the rows based on gene_id
        collected_data.drop_duplicates(subset=['gene_id'], inplace=True)

        # Populate gene_description and category using gene_id
        # Assuming gene_description is "gene_id tRNA"
        # Assuming category is simply "tRNA"
        collected_data['gene_description'] = collected_data['gene_id'].apply(lambda x: x + " tRNA")
        collected_data['category'] = "tRNA"
        collected_data['topic_ecosystem'] = "tRNA"
        collected_data['subcategory'] = ""  # Assuming subcategory is not defined

        # Count occurrences for each sample and fill the additional columns dynamically
        for sample in samples:
            collected_data[sample] = collected_data['gene_id'].map(lambda x: sample_counts[sample].count(x) if x in sample_counts[sample] else 0)

        # Sort the whole table by the gene_id column
        collected_data.sort_values(by='gene_id', inplace=True)

        # Write the collected data to the output file
        collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)
    """
}
