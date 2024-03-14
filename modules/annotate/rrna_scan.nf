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
            "--kingdom", "bacteria",
            fasta
        ]
        result = subprocess.run(barrnap_command, capture_output=True, text=True)
        if result.returncode != 0 or not result.stdout.strip():
            # Barrnap error or no output
            if verbose:
                print(f"Barrnap error or no rRNAs were detected for {sample_name}.", file=stderr)
            return None

        raw_rrna_table = pd.read_csv(
            io.StringIO(result.stdout),
            skiprows=1,
            sep="\t",
            header=None,
            names=RAW_RRNA_COLUMNS,
            index_col=0,
        )

        if raw_rrna_table.empty:
            # No rRNAs detected
            return None

        rrna_table_rows = []
        for gene, row in raw_rrna_table.iterrows():
            rrna_row_dict = {
                entry.split("=")[0]: entry.split("=")[1] for entry in row["note"].split(";")
            }
            rrna_table_rows.append([
                sample_name,
                row.begin,
                row.end,
                row.strand,
                rrna_row_dict["Name"].replace("_", " "),
                row["e-value"],
                rrna_row_dict.get("note", ""),
            ])

        return pd.DataFrame(
            rrna_table_rows,
            columns=RRNA_COLUMNS
        )

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

    rrna_df = run_barrnap("${fasta}", "${sample}", threads=${params.threads}, verbose=True)

    if rrna_df is not None and not rrna_df.empty:
        rrna_df.to_csv("${sample}_processed_rrnas.tsv", sep="\t", index=False)
    else:
        with open("${sample}_processed_rrnas.tsv", 'w') as output_file:
            output_file.write("NULL")
    """
}
