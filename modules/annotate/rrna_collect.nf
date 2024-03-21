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
    from collections import defaultdict
    import re

    def extract_rrna_gene_id(note):
        match = re.search(r'Name=(\\w+)_rRNA', note)
        if match:
            return f"{match.group(1)} rRNA"
        return 'Unknown rRNA'

    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]
    samples = [os.path.basename(f).replace('_processed_rrnas.tsv', '') for f in tsv_files]
    gene_type_counts = defaultdict(lambda: defaultdict(int))
    combined_data = []
    all_files_null = True

    for file in tsv_files:
        with open(file, 'r') as f:
            contents = f.read().strip()
            if contents == 'NULL':
                continue
            all_files_null = False
            df = pd.read_csv(file, sep='\\t')
            for index, row in df.iterrows():
                gene_id = extract_rrna_gene_id(row['note'])
                gene_type_counts[gene_id][row['sample']] += 1
            combined_data.append(df)

    if all_files_null:
        with open('collected_rrnas.tsv', 'w') as f: f.write('NULL')
        with open('combined_rrna_scan.tsv', 'w') as f: f.write('NULL')
    else:
        collected_data = []
        for gene_id, samples_counts in gene_type_counts.items():
            row = {'gene_id': gene_id, 'gene_description': f"{gene_id} gene", 'category': 'rRNA', 'topic_ecosystem': '', 'subcategory': ''}
            for sample in samples: row[sample] = samples_counts.get(sample, 0)
            collected_data.append(row)
        collected_df = pd.DataFrame(collected_data)
        collected_df.to_csv('collected_rrnas.tsv', sep='\\t', index=False)
        combined_df = pd.concat(combined_data, ignore_index=True)
        combined_df.to_csv('combined_rrna_scan.tsv', sep='\\t', index=False)
    """
}

