process TRNA_COLLECT {

    errorStrategy 'finish'

    input:
    file combined_trnas

    output:
    tuple val(sample), path("collected_trnas.tsv"), emit: trna_collected_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from ast import literal_eval

    def extract_samples_and_paths(combined_trnas):
        samples_and_paths = [(combined_trnas[i], combined_trnas[i + 1]) for i in range(0, len(combined_trnas), 2)]
        return samples_and_paths

    # Debugging statements
    print("Debug: Initial value of combined_trnas")
    print("${combined_trnas}")

    # Load and preprocess combined_trnas
    combined_trnas_list = literal_eval(open("${combined_trnas}").read())

    # Create an empty DataFrame to store the collected data
    collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "module", "header", "subheader"] + [sample for sample, _ in combined_trnas_list])

    # Iterate over each sample and corresponding path
    for sample, path in combined_trnas_list:
        try:
            with open(path, 'r') as file:
                # Read the processed tRNAs file for the current sample
                trna_data = pd.read_csv(file, sep="\t", skiprows=[0, 2])
                
                # Add data to the collected DataFrame
                # Update the following line based on how you want to populate the values
                # collected_data[sample] = ...

        except FileNotFoundError:
            print(f"Debug: File {path} not found.")
        except pd.errors.EmptyDataError:
            print(f"Debug: File {path} is empty.")
        except Exception as e:
            print(f"Debug: Error reading {sample}: {e}")

    # Write the collected data to the output file
    collected_data.to_csv("collected_trnas.tsv", sep="\t", index=False)

    """
}
