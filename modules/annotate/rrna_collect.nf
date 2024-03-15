process RRNA_COLLECT {

    errorStrategy 'finish'

    input:
    file combined_rrnas

    output:
    path("collected_rrnas.tsv"), emit: rrna_collected_out, optional: true
    path("combined_rrna_scan.tsv"), emit: rrna_combined_out, optional: true

    script:
    """
    import os
    import pandas as pd
    from collections import Counter

    # List all tsv files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Initialize a flag to check if all files are "NULL"
    all_files_null = True

    # Create a list to store individual input DataFrames
    individual_dfs = []
    # Create a list to store names of non-"NULL" files
    non_null_files = []

    for file in tsv_files:
        with open(file, 'r') as f:
            contents = f.read().strip()
            # Check if the file content is "NULL"
            if contents == "NULL":
                continue  # Skip this file
            all_files_null = False  # At least one file is not "NULL"
            input_df = pd.read_csv(file, sep='\t')
            individual_dfs.append(input_df)
            non_null_files.append(file)  # Add the file name to non_null_files

    # If all files contain "NULL", write "NULL" to output files and exit
    if all_files_null:
        with open("collected_rrnas.tsv", "w") as f:
            f.write("NULL")
        with open("combined_rrna_scan.tsv", "w") as f:
            f.write("NULL")
    else:
        # Extract sample names from non-null files
        samples = [os.path.basename(file).replace("_processed_rrnas.tsv", "") for file in non_null_files]
        collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "category", "topic_ecosystem", "subcategory"] + samples)
        sample_counts = {sample: Counter() for sample in samples}

        for df, file in zip(individual_dfs, non_null_files):  # Use non_null_files here
            sample_name = os.path.basename(file).replace("_processed_rrnas.tsv", "")
            gene_ids = df['type'].tolist()
            unique_gene_ids = list(set(gene_ids))
            gene_counts = Counter(gene_ids)

            collected_data = pd.concat([collected_data, pd.DataFrame({'gene_id': unique_gene_ids})], ignore_index=True)
            collected_data['category'] = 'rRNA'
            for gene_id, count in gene_counts.items():
                sample_counts[sample_name][gene_id] += count

        combined_df = pd.concat(individual_dfs, ignore_index=True)
        combined_df.to_csv("combined_rrna_scan.tsv", sep="\t", index=False)

        collected_data['type'] = ""
        collected_data['gene_description'] = collected_data['gene_id'] + " gene"
        collected_data['topic_ecosystem'] = ""
        collected_data['subcategory'] = ""
        collected_data.drop_duplicates(subset=['gene_id'], inplace=True)

        for sample in samples:
            collected_data[sample] = collected_data['gene_id'].map(lambda x: sample_counts[sample][x] if x in sample_counts[sample] else 0)

        collected_data = collected_data[collected_data['gene_id'] != 'gene_id']
        collected_data.sort_values(by='gene_id', inplace=True)
        collected_data.drop(['type'], axis=1, inplace=True)
        collected_data.to_csv("collected_rrnas.tsv", sep="\t", index=False)

    """
}
