process TRNA_COLLECT {

    errorStrategy 'finish'

    input:
    val combined_trnas

    output:
    tuple val(sample), path("collected_trnas.tsv"), emit: trna_collected_out, optional: true

    script:
    """
    #!/usr/bin/env python

    # Replace single quotes with double quotes in the combined_trnas variable
    combined_trnas = """${combined_trnas.replaceAll("'", '\"')}"""

    import pandas as pd

    # Function to read the first two lines of each input file
    def read_first_two_lines(file_path):
        with open(file_path, 'r') as file:
            lines = [file.readline().strip() for _ in range(2)]
        return lines

    # Split the modified string into a list
    combined_trnas_list = eval(combined_trnas)

    # Iterate over samples and file paths in combined_trnas
    for i in range(0, len(combined_trnas_list), 2):
        sample = combined_trnas_list[i]
        file_path = combined_trnas_list[i + 1]

        # Read the first two lines of the input file
        lines = read_first_two_lines(file_path)

        # Print sample name, file path, and first two lines
        print(f"Sample: {sample}, File: {file_path}")
        print("\n".join(lines))
        print("\n" + "=" * 50 + "\n")  # Separator line
    """
}
