process RRNA_SCAN {
    label 'process_low'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_barrnap_trnascan-se:ed2ab26abf39304b"

    tag { input_fasta }

    input:
    tuple val(input_fasta), path(fasta)

    output:
    tuple val(input_fasta), path("${input_fasta}_processed_rrnas.tsv"), emit: rrna_scan_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import subprocess
    import io
    from sys import stderr

    def run_barrnap(fasta, input_fasta_name, threads, verbose=True):
        barrnap_command = [
            "barrnap",
            "--threads", str(threads),
            "--kingdom", "bacteria",  # Adjust as necessary
            fasta
        ]
        result = subprocess.run(barrnap_command, capture_output=True, text=True)
        raw_rrna_str = result.stdout

        if not raw_rrna_str.strip():
            print(f"No rRNAs were detected for {input_fasta_name}.", file=stderr)
            return pd.DataFrame(columns=["input_fasta", "query_id", "type", "begin", "end", "strand", "e-value", "note"])  # Ensure this matches RRNA_COLUMNS

        try:
            rrna_df = pd.read_csv(
                io.StringIO(raw_rrna_str),
                sep="\t",
                header=None,
                names=["query_id", "tool_name", "type", "begin", "end", "strand", "e-value", "score", "note"],
                usecols=["query_id", "type", "begin", "end", "strand", "e-value", "note"],
                comment='#'  # This will skip lines starting with '#', including the '##gff-version 3' line
            )
            rrna_df.insert(0, 'input_fasta', input_fasta_name)
            return rrna_df
        except pd.errors.ParserError:
            print(f"Parser error processing barrnap output for {input_fasta_name}. Output may not be in the expected format.", file=stderr)
            return pd.DataFrame(columns=["input_fasta", "query_id", "type", "begin", "end", "strand", "e-value", "note"])


    rrna_df = run_barrnap("${fasta}", "${input_fasta}", threads=${params.threads}, verbose=True)

    if not rrna_df.empty:
        rrna_df.to_csv(f"${input_fasta}_processed_rrnas.tsv", sep="\t", index=False)
    else:
        with open(f"${input_fasta}_processed_rrnas.tsv", "w") as file:
            file.write("NULL")

    """
}
