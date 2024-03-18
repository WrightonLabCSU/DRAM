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
    from collections import Counter, defaultdict

    # List all TSV files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Extract sample names from the file names
    samples = [os.path.basename(file).replace("_processed_rrnas.tsv", "") for file in tsv_files]

    # Check if all files contain "NULL"
    all_files_null = all(open(file).read().strip() == "NULL" for file in tsv_files)

    if all_files_null:
        # Write "NULL" to output files and exit
        with open("collected_rrnas.tsv", "w") as f:
            f.write("NULL")
        with open("combined_rrna_scan.tsv", "w") as f:
            f.write("NULL")
    else:
        # Initialize a dictionary to store counts of each specific rRNA gene type per sample
        gene_type_counts = defaultdict(lambda: defaultdict(int))

        # Initialize list for storing data from non-"NULL" files for combined output
        combined_data = []

        for file in tsv_files:
            with open(file, 'r') as f:
                contents = f.read().strip()
                if contents == "NULL":
                    continue  # Skip this file

                df = pd.read_csv(file, sep='\t')
                # Update combined_data
                combined_data.extend(df.to_dict('records'))

                for _, row in df.iterrows():
                    gene_type = row['type']  # Extract specific rRNA gene type
                    gene_type_counts[gene_type][row['sample']] += 1

        # Prepare collected data DataFrame from gene_type_counts
        collected_data_rows = []
        for gene_type, samples_counts in gene_type_counts.items():
            row = {
                'gene_id': gene_type,
                'gene_description': f"{gene_type} gene",
                'category': 'rRNA',
                'topic_ecosystem': '',
                'subcategory': ''
            }
            # Update counts for each sample
            for sample in samples:
                row[sample] = samples_counts.get(sample, 0)
            collected_data_rows.append(row)

        collected_df = pd.DataFrame(collected_data_rows)
        collected_df = collected_df[['gene_id', 'gene_description', 'category', 'topic_ecosystem', 'subcategory'] + samples]
        collected_df.sort_values(by='gene_id', inplace=True)
        collected_df.to_csv("collected_rrnas.tsv", sep='\t', index=False)

        # Combine non-"NULL" data for combined_rrna_scan.tsv
        combined_df = pd.DataFrame(combined_data)
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