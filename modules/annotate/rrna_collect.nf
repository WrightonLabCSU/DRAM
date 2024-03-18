process RRNA_COLLECT {

    errorStrategy 'finish'

    input:
    file combined_rrnas

    output:
    path("collected_rrnas.tsv"), emit: rrna_collected_out, optional: true
    path("combined_rrna_scan.tsv"), emit: rrna_combined_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import os
    import pandas as pd
    from collections import Counter

    # List all TSV files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Extract sample names from the file names
    samples = [os.path.basename(file).replace("_processed_rrnas.tsv", "") for file in tsv_files]

    # Check if all files contain "NULL"
    all_files_null = all([open(file).read().strip() == "NULL" for file in tsv_files])

    if all_files_null:
        # Write "NULL" to output files and exit
        with open("collected_rrnas.tsv", "w") as f:
            f.write("NULL")
        with open("combined_rrna_scan.tsv", "w") as f:
            f.write("NULL")
    else:
        # Initialize a DataFrame to store aggregated data
        collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "category", "topic_ecosystem", "subcategory"] + samples)

        # Initialize a Counter for each sample to track gene counts
        sample_counts = {sample: Counter() for sample in samples}

        # Initialize a list to store data from non-"NULL" files
        individual_dfs = []

        # Iterate through each input file
        for file in tsv_files:
            with open(file, 'r') as f:
                contents = f.read().strip()
                if contents == "NULL":
                    continue  # Skip this file
                
                sample_name = os.path.basename(file).replace("_processed_rrnas.tsv", "")
                df = pd.read_csv(file, sep='\t')

                # Increment count for each gene_id in this sample
                for gene_id in df['type'].unique():
                    sample_counts[sample_name][gene_id] += len(df[df['type'] == gene_id])
                
                individual_dfs.append(df)

        # Populate collected_data based on aggregated sample_counts
        for sample, counts in sample_counts.items():
            for gene_id, count in counts.items():
                if gene_id not in collected_data['gene_id'].values:
                    collected_data = collected_data.append({
                        'gene_id': gene_id,
                        'gene_description': f"{gene_id} gene",
                        'category': 'rRNA',
                        'topic_ecosystem': '',
                        'subcategory': ''
                    }, ignore_index=True)
                collected_data.loc[collected_data['gene_id'] == gene_id, sample] = count

        # Ensure correct order of columns and sort by 'gene_id'
        collected_data = collected_data[['gene_id', 'gene_description', 'category', 'topic_ecosystem', 'subcategory'] + samples]
        collected_data.sort_values(by='gene_id', inplace=True)
        collected_data.to_csv("collected_rrnas.tsv", sep='\t', index=False)

        # Combine non-"NULL" dataframes for combined_rrna_scan.tsv, ensuring 'query_id' is included
        combined_df = pd.concat(individual_dfs, ignore_index=True)
        combined_df.to_csv("combined_rrna_scan.tsv", sep='\t', index=False)

    """
}


/* OLD VERSION
process RRNA_COLLECT {

    errorStrategy 'finish'

    input:
    file combined_rrnas

    output:
    path("collected_rrnas.tsv"), emit: rrna_collected_out, optional: true
    path("combined_rrna_scan.tsv"), emit: rrna_combined_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import os
    import pandas as pd
    from collections import Counter

    # List all tsv files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Extract sample names from the file names
    samples = [os.path.basename(file).replace("_processed_rrnas.tsv", "") for file in tsv_files]

    # Create an empty DataFrame to store the collected data
    collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "category", "topic_ecosystem", "subcategory"] + samples)

    # Create a dictionary to store counts for each sample
    sample_counts = {sample: Counter() for sample in samples}

    # Create a list to store individual input DataFrames
    individual_dfs = []

    # Iterate through each input file
    for file in tsv_files:
        sample_name = os.path.basename(file).replace("_processed_rrnas.tsv", "")
        input_df = pd.read_csv(file, sep='\t')

        # Populate gene_id column with collective values from the "type" column
        gene_ids = input_df['type'].tolist()
        unique_gene_ids = list(set(gene_ids))
        gene_counts = Counter(gene_ids)

        # Populate gene_description column
        collected_data = pd.concat([collected_data, pd.DataFrame({'gene_id': unique_gene_ids})], ignore_index=True)

        # Set category column values to "rRNA"
        collected_data['category'] = 'rRNA'

        # Update counts for each sample
        for gene_id, count in gene_counts.items():
            sample_counts[sample_name][gene_id] += count

        # Append the individual input DataFrame to the list
        individual_dfs.append(input_df)

    # Combine all individual DataFrames into a single DataFrame
    combined_df = pd.concat(individual_dfs, ignore_index=True)

    # Write the combined DataFrame to the output file
    combined_df.to_csv("combined_rrna_scan.tsv", sep="\t", index=False)

    # Add 'type' column to collected_data
    collected_data['type'] = ""

    # Populate other columns based on the given rules
    collected_data['gene_description'] = collected_data['gene_id'] + " gene"
    collected_data['topic_ecosystem'] = ""
    collected_data['subcategory'] = ""

    # Deduplicate the rows based on gene_id
    collected_data.drop_duplicates(subset=['gene_id'], inplace=True)

    # Count occurrences for each sample and fill the additional columns dynamically
    for sample in samples:
        collected_data[sample] = collected_data['gene_id'].map(lambda x: sample_counts[sample][x] if x in sample_counts[sample] else 0)

    # Remove the first row that includes extra column names
    collected_data = collected_data[collected_data['gene_id'] != 'gene_id']

    # Sort the whole table by the gene_id column
    collected_data.sort_values(by='gene_id', inplace=True)

    # Drop the 'type' column
    collected_data.drop(['type'], axis=1, inplace=True)

    # Write the collected data to the output file
    collected_data.to_csv("collected_rrnas.tsv", sep="\t", index=False)


    """
}
*/