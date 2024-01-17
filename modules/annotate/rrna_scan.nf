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

    def run_barrnap(fasta, fasta_name, logger, threads=10, verbose=True):
        # The provided run_process function is assumed to be defined elsewhere in your code
        raw_rrna_str = run_process(
            ["barrnap", "--threads", str(threads), fasta],
            logger,
            capture_stdout=True,
            check=False,
            verbose=verbose,
        )
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
                    fasta_name,
                    row.begin,
                    row.end,
                    row.strand,
                    rrna_row_dict["Name"].replace("_", " "),
                    row["e-value"],
                    rrna_row_dict.get("note", ""),
                ]
            )
        if len(raw_rrna_table) > 0:
            return pd.DataFrame(
                rrna_table_rows, index=raw_rrna_table.index, columns=RRNA_COLUMNS
            ).reset_index()
        else:
            logger.warning("No rRNAs were detected, no rrnas.tsv file will be created.")
            return None

    RAW_RRNA_COLUMNS = [
        "scaffold",
        "tool_name",
        "type",
        "begin",
        "end",
        "e-value",
        "strand",
        "empty",
        "note",
    ]
    RRNA_COLUMNS = ["fasta", "begin", "end", "strand", "type", "e-value", "note"]

    # Run barrnap
    rrna_df = run_barrnap("${fasta}", "${sample}", logger, threads=${params.threads}, verbose=True)

    if rrna_df is not None:
        # Save rrna_df to ${sample}_processed_rrnas.tsv
        rrna_df.to_csv("${sample}_processed_rrnas.tsv", sep="\t", index=False)
    else:
        # No rRNAs detected, create an empty file
        with open("${sample}_processed_rrnas.tsv", 'w') as empty_file:
            pass
    """
}
