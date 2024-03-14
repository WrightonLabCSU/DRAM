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

    def check_file_content(file_path):
        with open(file_path, 'r') as file:
            content = file.read().strip()
            return content == "NULL"

    # List all tsv files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Check if file contents are "NULL"
    all_files_null = all(check_file_content(f) for f in tsv_files)

    if all_files_null:
        # Write "NULL" to both output files if all files contain "NULL"
        with open("combined_rrna_scan.tsv", "w") as combined_out, open("collected_rrnas.tsv", "w") as collected_out:
            combined_out.write("NULL")
            collected_out.write("NULL")
    else:
        # Process files normally, skipping any files that contain "NULL"
        tsv_files = [f for f in tsv_files if not check_file_content(f)]

        samples = [os.path.basename(file).replace("_processed_rrnas.tsv", "") for file in tsv_files]

        collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "category", "topic_ecosystem", "subcategory"] + samples)

        sample_counts = {sample: Counter() for sample in samples}

        individual_dfs = []

        for file in tsv_files:
            sample_name = os.path.basename(file).replace("_processed_rrnas.tsv", "")
            input_df = pd.read_csv(file, sep='\t')

            gene_ids = input_df['type'].tolist()
            unique_gene_ids = list(set(gene_ids))
            gene_counts = Counter(gene_ids)

            collected_data = pd.concat([collected_data, pd.DataFrame({'gene_id': unique_gene_ids})], ignore_index=True)

            collected_data['category'] = 'rRNA'

            for gene_id, count in gene_counts.items():
                sample_counts[sample_name][gene_id] += count

            individual_dfs.append(input_df)

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
