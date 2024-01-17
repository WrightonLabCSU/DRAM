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

    def extract_samples_and_paths(combined_trnas):
        samples_and_paths = []
        for i in range(0, len(combined_trnas), 2):
            sample = combined_trnas[i].strip('[], ')
            path = combined_trnas[i + 1].strip('[], ')
            samples_and_paths.append((sample, path))
        return samples_and_paths

    # Debugging statements
    print("Debug: Before extracting samples and paths")
    combined_trnas_str = "input.1 bin-124_processed_trnas.tsv input.2 bin-129_processed_trnas.tsv"
    print(f"Debug: Combined_trnas_str: {combined_trnas_str}")
    combined_trnas_list = extract_samples_and_paths(combined_trnas_str)
    print(f"Debug: Combined_trnas_list: {combined_trnas_list}")

    # Print the first two lines of each input file
    for sample, path in combined_trnas_list:
        try:
            with open(path, 'r') as file:
                lines = file.readlines()[:2]
                print(f"Debug: First two lines of {sample}: {lines}")
        except FileNotFoundError:
            print(f"Debug: File {path} not found.")
        except Exception as e:
            print(f"Debug: Error reading {sample}: {e}")


    """
}
