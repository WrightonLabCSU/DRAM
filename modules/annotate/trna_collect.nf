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

    # List all tsv files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Check if all files contain "NULL"
    all_files_null = all(open(file).read().strip() == "NULL" for file in tsv_files)

    if all_files_null:
        # If all files contain "NULL", write "NULL" to the output and exit
        with open("collected_trnas.tsv", "w") as output_file:
            output_file.write("NULL")
    else:
        # Process files that do not contain "NULL"
        samples = [os.path.basename(file).replace("_processed_trnas.tsv", "") for file in tsv_files if open(file).read().strip() != "NULL"]

        # Create an empty DataFrame to store the collected data
        collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "category", "topic_ecosystem", "subcategory"] + samples)

        # Create a dictionary to store counts for each sample
        sample_counts = {sample: [] for sample in samples}

        # Iterate through each input file, skipping files with "NULL" content
        for file in tsv_files:
            if open(file).read().strip() == "NULL":
                continue  # Skip this file
            
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
        collected_data['gene_description'] = collected_data['gene_id'] + " tRNA with " + collected_data['codon'] + " Codon"
        collected_data['category'] = collected_data['gene_id'] + " tRNA"

        # Remove parentheses and text in between from gene_description and category columns
        collected_data['gene_description'] = collected_data['gene_id'].apply(lambda x: re.split(r'\(|\)', x)[0].strip() + " tRNA with " + collected_data.loc[collected_data['gene_id'] == x, 'codon'].values[0] + " Codon")
        collected_data['category'] = collected_data['gene_id'].apply(lambda x: re.split(r'\(|\)', x)[0].strip() + " tRNA")

        collected_data['topic_ecosystem'] = "tRNA"
        collected_data['subcategory'] = ""

        # Deduplicate the rows based on gene_id
        collected_data.drop_duplicates(subset=['gene_id'], inplace=True)

        # Count occurrences for each sample and fill the additional columns dynamically
        for sample in samples:
            collected_data[sample] = collected_data['gene_id'].map(lambda x: sample_counts[sample].count(x) if x in sample_counts[sample] else 0)

        # Remove 'type' and 'codon' columns
        collected_data.drop(['type', 'codon'], axis=1, inplace=True)

        # Remove the first row that includes extra column names
        collected_data = collected_data[collected_data['gene_id'] != 'gene_id']

        # Sort the whole table by the gene_id column
        collected_data.sort_values(by='gene_id', inplace=True)

        # Write the collected data to the output file
        collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)

    """
}
