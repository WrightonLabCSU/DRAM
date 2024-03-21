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

    print("Starting TRNA_COLLECT script...")

    # List all tsv files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]
    print(f"Found {len(tsv_files)} TSV files.")

    # Initialize variable to check if all files are "NULL"
    all_files_null = True

    # Check each file
    for file in tsv_files:
        with open(file, 'r') as f:
            contents = f.read().strip()
            if contents != "NULL":
                all_files_null = False
                break

    if all_files_null:
        print("All files contain 'NULL'. Writing 'NULL' to collected_trnas.tsv")
        with open('collected_trnas.tsv', 'w') as f:
            f.write("NULL")
    else:
        print("Processing TSV files...")

        # Initialize an empty DataFrame for collected data
        collected_data = pd.DataFrame()

        # Process each file
        for file in tsv_files:
            print(f"Processing file: {file}")
            # Skip files with 'NULL' content
            with open(file, 'r') as f:
                if f.read().strip() == "NULL":
                    continue
            
            # Read the TSV file
            df = pd.read_csv(file, sep='\t')
            # Assuming the structure of your TSV files matches the example given
            # Append to the collected_data DataFrame
            collected_data = pd.concat([collected_data, df], ignore_index=True)

        # Assuming you want to perform some aggregations or transformations
        # For demonstration, I'll simply remove duplicates based on 'gene_id'
        if not collected_data.empty:
            collected_data = collected_data.drop_duplicates(subset=['gene_id'])

            # Write the processed data to a TSV file
            collected_data.to_csv('collected_trnas.tsv', sep='\t', index=False)
            print("Processed data written to collected_trnas.tsv")
        else:
            print("No data to process after checking files.")

    """
}
