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
        if not raw_rrna_str.strip():
            # No output from barrnap, indicating no rRNAs were detected
            if verbose:
                print(f"No rRNAs were detected for {sample_name}.", file=stderr)
            return None
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
            rrna_table_rows, index=raw_rrna_table.index, columns=RRNA_COLUMNS
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

    if rrna_df is not None:
        # Save rrna_df to ${sample}_processed_rrnas.tsv
        rrna_df.to_csv("${sample}_processed_rrnas.tsv", sep="\t", index=False)
    else:
        # Copy the distill_dummy_sheet content to the output file
        # Assumes that distill_dummy_sheet is a path to the file with default contents
        with open(${params.distill_dummy_sheet}, 'r') as dummy_sheet:
            with open("${sample}_processed_rrnas.tsv", 'w') as output_file:
                output_file.write(dummy_sheet.read())
    """
}
