process RRNA_SCAN {

    errorStrategy 'finish'

    tag { sample }

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_processed_rrnas.tsv"), emit: rrna_scan_out, optional: true

    script:

    """
    RAW_RRNA_COLUMNS=(
        "scaffold"
        "tool_name"
        "type"
        "begin"
        "end"
        "e-value"
        "strand"
        "empty"
        "note"
    )
    
    RRNA_COLUMNS=(
        "fasta"
        "begin"
        "end"
        "strand"
        "type"
        "e-value"
        "note"
    )

    raw_rrna_str=\$(barrnap --threads ${params.threads} ${fasta})
    
    # Process barrnap Output
    if [ -n "\$raw_rrna_str" ]; then
        raw_rrna_table=\$(echo "\$raw_rrna_str" | tail -n +2 | awk -F'\t' '{print \$1, \$2, \$3, \$4, \$5, \$6, \$7}')

        rrna_table_rows=()
        while read -r \${RAW_RRNA_COLUMNS[@]}; do
            rrna_row_dict=\$(awk -F';' '{ for (i=1; i<=NF; i++) print \$i }' <<< "\$note" | awk -F'=' '{print \$1, \$2}')

            rrna_table_rows+=("\${sample}\t\$begin\t\$end\t\$strand\t\${rrna_row_dict["Name"]//_/ }\t\$evalue\t\${rrna_row_dict["note"]}")
        done <<< "\$raw_rrna_table"

        # Create DataFrame from the processed rRNA table rows
        rrna_frame=\$(awk 'BEGIN{OFS="\t"}{print \$0}' <<< "\${rrna_table_rows[@]}")

        # Save the processed rRNA frame to a new file
        echo -e "\$rrna_frame" > ${sample}_processed_rrnas.tsv
    else
        echo "No rRNAs were detected, no rrnas.tsv file will be created."
        exit 1
    fi
    """

}
