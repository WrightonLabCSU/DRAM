process RRNA_SCAN {

    errorStrategy 'finish'

    tag { sample }

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_processed_rrnas.tsv"), emit: rrna_scan_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import subprocess
    import io
    from sys import stderr

    def run_barrnap(fasta, sample_name, threads, verbose=True):
        barrnap_command = [
            "barrnap",
            "--threads", str(threads),
            "--kingdom", "bacteria",  # Add any other necessary barrnap options
            fasta
        ]
        raw_rrna_str = subprocess.run(barrnap_command, capture_output=True, text=True).stdout
        if not raw_rrna_str.strip():  # Check if barrnap output is empty
            print(f"No rRNAs were detected for {sample_name}.", file=stderr)
            return pd.DataFrame(columns=["sample", "query_id", "type", "begin", "end", "strand", "e-value", "note"])  # Ensure this matches RRNA_COLUMNS

        # Attempt to read the output into a DataFrame, assuming the output includes headers or using provided column names
        try:
            rrna_df = pd.read_csv(
                io.StringIO(raw_rrna_str),
                sep="\t",
                header=None,
                names=["query_id", "tool_name", "type", "begin", "end", "strand", "e-value", "score", "note"],
                usecols=["query_id", "type", "begin", "end", "strand", "e-value", "note"]
            )
            rrna_df.insert(0, 'sample', sample_name)  # Add 'sample' column at the beginning
            return rrna_df
        except pd.errors.EmptyDataError:
            print(f"Error processing barrnap output for {sample_name}.", file=stderr)
            return pd.DataFrame(columns=["sample", "query_id", "type", "begin", "end", "strand", "e-value", "note"])  # Match to expected columns

    # Using parameters from Nextflow
    fasta = "${fasta}"
    sample_name = "${sample}"
    threads = ${params.threads}

    # Run barrnap
    rrna_df = run_barrnap(fasta, sample_name, threads, verbose=True)

    if not rrna_df.empty:
        rrna_df.to_csv(f"{sample_name}_processed_rrnas.tsv", sep="\t", index=False)
    else:
        # Write "NULL" to the file if rrna_df is empty or in case of an error
        with open(f"{sample_name}_processed_rrnas.tsv", "w") as file:
            file.write("NULL")

    """
}
