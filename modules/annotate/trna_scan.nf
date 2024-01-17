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

    def process_trnascan_output(input_file, output_file, sample_name):
        # Read the input file into a DataFrame
        trna_frame = pd.read_csv(input_file, sep="\t", skiprows=[0, 2])

        # Strip leading and trailing spaces from column names
        trna_frame.columns = trna_frame.columns.str.strip()

        # Add a new "sample" column and populate it with the sample_name value
        trna_frame.insert(0, "sample", sample_name)

        # Remove the "Note" column if present
        trna_frame = trna_frame.drop(columns=["Note"], errors="ignore")

        # Keep only the first occurrence of "Begin" and "End" columns
        trna_frame = trna_frame.loc[:, ~trna_frame.columns.duplicated(keep='first')]

        # Reorder columns
        columns_order = ["sample", "Name", "tRNA #", "Begin", "End", "Type", "Codon", "Score"]

        # Rename specified columns
        trna_frame = trna_frame.rename(columns={"Name": "query_id", "Begin": "begin", "End": "end", "Type": "type", "Codon": "codon", "Score": "score"})

        # Write the processed DataFrame to the output file
        trna_frame.to_csv(output_file, sep="\t", index=False)

    # Run tRNAscan-SE
    trna_out = "${sample}_trna_out.txt"
    subprocess.run(["tRNAscan-SE", "-G", "-o", trna_out, "--thread", "${params.threads}", "${fasta}"], check=True)

    # Process tRNAscan-SE output
    process_trnascan_output(trna_out, "${sample}_processed_trnas.tsv", "${sample}")
    """
}
