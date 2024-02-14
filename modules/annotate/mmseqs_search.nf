process MMSEQS_SEARCH {

    tag { sample }

    input:
    tuple val( sample ), path( query_database, stageAs: "query_database/" )
    path( mmseqs_database )
    val( bit_score_threshold)
    file( db_descriptions )
    val( db_name )

    output:
    tuple val( sample ), path("mmseqs_out/${sample}_mmseqs_${db_name}.tsv"), emit: mmseqs_search_out
    tuple val( sample ), path("mmseqs_out/${sample}_mmseqs_${db_name}_formatted.csv"), emit: mmseqs_search_formatted_out

    script:
    """
    #!/usr/bin/env python
    import os
    import subprocess
    import pandas as pd

    # Create symbolic link to mmseqs database
    os.symlink(str(mmseqs_database), './')

    # Create temporary directory
    os.makedirs('mmseqs_out/tmp', exist_ok=True)

    print(f"Query Database: {query_database}")
    print(f"Target Database: {mmseqs_database}")

    # Perform search
    subprocess.run([
        'mmseqs', 'search', f'query_database/{sample}.mmsdb',
        f'{db_name}.mmsdb', f'mmseqs_out/{sample}_{db_name}.mmsdb',
        'mmseqs_out/tmp', '--threads', str(params.threads)
    ])

    # Filter to only best hit
    subprocess.run([
        'mmseqs', 'filterdb', f'mmseqs_out/{sample}_{db_name}.mmsdb',
        f'mmseqs_out/{sample}_{db_name}_tophit.mmsdb', '--extract-lines', '1'
    ])

    # Filter to only hits with minimum bit score
    subprocess.run([
        'mmseqs', 'filterdb', '--filter-column', '2', '--comparison-operator', 'ge',
        '--comparison-value', str(bit_score_threshold), '--threads', str(params.threads),
        f'mmseqs_out/{sample}_{db_name}_tophit.mmsdb',
        f'mmseqs_out/{sample}_{db_name}_tophit_minbitscore{bit_score_threshold}.mmsdb'
    ])

    # Convert results to BLAST outformat 6
    subprocess.run([
        'mmseqs', 'convertalis', f'query_database/{sample}.mmsdb',
        f'{db_name}.mmsdb', f'mmseqs_out/{sample}_{db_name}_tophit_minbitscore{bit_score_threshold}.mmsdb',
        f'mmseqs_out/{sample}_mmseqs_{db_name}.tsv', '--threads', str(params.threads)
    ])

    # Define the input and output file paths
    input_path = f"mmseqs_out/{sample}_mmseqs_{db_name}.tsv"
    output_path = f"mmseqs_out/{sample}_mmseqs_{db_name}_formatted.csv"

    # Process MMseqs output
    df_mmseqs = pd.read_csv(input_path, sep='\t', header=None, names=['query_id', 'target_id', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    df_mmseqs['start_position'] = df_mmseqs['qstart']
    df_mmseqs['end_position'] = df_mmseqs['qend']
    df_mmseqs[f"{db_name}_id"] = df_mmseqs['target_id']
    df_mmseqs[f"{db_name}_bitScore"] = df_mmseqs['bitscore']
    df_mmseqs_final = df_mmseqs[['query_id', 'start_position', 'end_position', f"{db_name}_id", f"{db_name}_bitScore"]]

    # Check if db_descriptions contains 'NULL'
    with open(str(db_descriptions), 'r') as file:
        if file.read().strip() == 'NULL':
            df_mmseqs_final.to_csv(output_path, index=False)
        else:
            # Reload db_descriptions to process it
            file.seek(0)  # Reset file pointer to the beginning
            df_descriptions = pd.read_csv(str(db_descriptions), sep='\t', header=None)
            # Dynamically generate column names based on db_name and existing column indexes
            desc_columns = [f"{db_name}_{i}" for i in range(1, len(df_descriptions.columns))]
            df_descriptions.columns = ['target_id'] + desc_columns
            
            # Merge MMseqs output with descriptions based on the target_id
            df_merged = pd.merge(df_mmseqs_final, df_descriptions, left_on=f"{db_name}_id", right_on='target_id', how='left')
            # Drop the redundant 'target_id' column from the merge if needed
            df_merged.drop(columns=['target_id'], inplace=True)
            df_merged.to_csv(output_path, index=False)

    """
}
