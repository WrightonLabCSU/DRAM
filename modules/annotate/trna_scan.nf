process TRNA_SCAN {

    errorStrategy 'finish'

    tag { sample }

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_processed_trnas.tsv"), emit: trna_scan_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import subprocess

    def process_trnascan_output(input_file, output_file):
        # Read the input file into a DataFrame
        trna_frame = pd.read_csv(input_file, sep="\t", skiprows=[0, 2])

        # Print column names for debugging
        print("Column names before processing:")
        print(trna_frame.columns)

        # Remove the "Note" column if present
        trna_frame = trna_frame.drop(columns=["Note"], errors="ignore")

        # Remove the first and third lines
        trna_frame = trna_frame.iloc[1::2]

        # Keep only the first occurrence of "Begin" and "End" columns
        unique_columns = trna_frame.columns.difference(["Begin", "End"])
        trna_frame = trna_frame.loc[:, unique_columns].join(trna_frame[["Begin", "End"]].apply(lambda col: col.first_valid_index(), axis=1).rename(columns={0: "BeginEndIndex"}))

        # Print column names again for debugging
        print("Column names after processing:")
        print(trna_frame.columns)

        # Reorder columns
        columns_order = ["Name", "tRNA #", "Begin", "End", "Type", "Codon", "Score"]

        # Print debugging information
        print("Columns to keep in order:")
        print(columns_order)

        # Try to reorder columns
        try:
            trna_frame = trna_frame[columns_order]
        except KeyError as e:
            print(f"Error: {e}")

        # Write the processed DataFrame to the output file
        trna_frame.to_csv(output_file, sep="\t", index=False)

    # Run tRNAscan-SE
    trna_out = "${sample}_trna_out.txt"
    subprocess.run(["tRNAscan-SE", "-G", "-o", trna_out, "--thread", "${params.threads}", "${fasta}"], check=True)

    # Process tRNAscan-SE output
    process_trnascan_output(trna_out, "${sample}_processed_trnas.tsv")
    """
}
