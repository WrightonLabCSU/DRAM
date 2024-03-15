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
    import io
    import subprocess
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
            return pd.DataFrame()  # Return an empty DataFrame

        raw_rrna_table = pd.read_csv(
            io.StringIO(raw_rrna_str),
            skiprows=1,
            sep="\t",
            header=None,
            names=RAW_RRNA_COLUMNS,
            index_col=0,
        )
        rrna_table_rows = list()
        for gene, row in raw_rrna_table.iterrows():
            rrna_row_dict = {
                entry.split("=")[0]: entry.split("=")[1] for entry in row["note"].split(";")
            }
            rrna_table_rows.append(
                [
                    sample_name,
                    row.begin,
                    row.end,
                    row.strand,
                    rrna_row_dict["Name"].replace("_", " "),
                    row["e-value"],
                    rrna_row_dict.get("note", ""),
                ]
            )
        return pd.DataFrame(
            rrna_table_rows, columns=RRNA_COLUMNS
        ).reset_index(drop=True)

    RAW_RRNA_COLUMNS = [
        "query_id",
        "tool_name",
        "type",
        "begin",
        "end",
        "e-value",
        "strand",
        "empty",
        "note",
    ]
    RRNA_COLUMNS = ["sample", "begin", "end", "strand", "type", "e-value", "note"]

    # Run barrnap
    rrna_df = run_barrnap("${fasta}", "${sample}", threads=${params.threads}, verbose=True)

    if rrna_df is not None and not rrna_df.empty:
        rrna_df.to_csv("${sample}_processed_rrnas.tsv", sep="\t", index=False)
    else:
        # Write "NULL" to the file if rrna_df is empty
        with open("${sample}_processed_rrnas.tsv", "w") as file:
            file.write("NULL")

    """
}
