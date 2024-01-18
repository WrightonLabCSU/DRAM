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

    # Extract sample names from the file names
    samples = [os.path.basename(file).replace("_processed_trnas.tsv", "") for file in tsv_files]

    # Create an empty DataFrame to store the collected data
    collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "module", "header", "subheader"] + samples)

    # Iterate over each input file
    for file, sample in zip(tsv_files, samples):
        try:
            print(f"Debug: Processing {file} for sample {sample}")

            # Read the processed tRNAs file for the current sample
            trna_data = pd.read_csv(file, sep="\t")

            # Debugging statement
            print(f"Debug: Data for {sample}:\n{trna_data}")

            # Generate values for each column based on the clarified requirements
            gene_id = f"{trna_data['type'].iloc[0]} ({trna_data['codon'].iloc[0]})"
            gene_description = f"{trna_data['type'].iloc[0]} tRNA with {trna_data['codon'].iloc[0]} Codon"
            module = f"{trna_data['type'].iloc[0]} tRNA"
            header = "tRNA"
            subheader = ""

            # Debugging statements
            print(f"Debug: gene_id for {sample}: {gene_id}")
            print(f"Debug: gene_description for {sample}: {gene_description}")
            print(f"Debug: module for {sample}: {module}")
            print(f"Debug: header for {sample}: {header}")
            print(f"Debug: subheader for {sample}: {subheader}")

            # Add data to the collected DataFrame
            collected_data[sample] = 0

        except FileNotFoundError:
            print(f"Debug: File {file} not found.")
        except pd.errors.EmptyDataError:
            print(f"Debug: File {file} is empty.")
        except Exception as e:
            print(f"Debug: Error reading {sample}: {e}")

    # Debugging statement
    print(f"Debug: Final collected_data DataFrame:\n{collected_data}")

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)

    """
}
