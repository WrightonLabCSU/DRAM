process TRNA_SCAN {

    errorStrategy 'finish'

    tag { sample }

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_processed_trnas.tsv"), emit: trna_scan_out, optional: true

    script:

    script:

    """
    tRNAscan-SE \\
    -G \\
    --thread ${params.threads} \\
    -o ${sample}_trna_out.txt \\
    ${fasta}

    # Process tRNAscan-SE Output

    # Remove first and third lines
    sed -i '1d;3d' ${sample}_trna_out.txt

    # Remove "Note" column
    awk '{NF=NF-1}1' ${sample}_trna_out.txt > ${sample}_trna_out_temp.txt
    mv ${sample}_trna_out_temp.txt ${sample}_trna_out.txt

    # Remove extra occurrences of "Begin" and "End" columns
    awk '!seen[\$3,\$4]++ {print \$1, \$2, \$5, \$6, \$7, \$8, \$9}' ${sample}_trna_out.txt > ${sample}_processed_trnas_temp.tsv

    # Extract the "Score" column
    awk '!seen[\$3,\$4,\$9]++ {print \$9}' ${sample}_trna_out.txt > ${sample}_score_column.txt

    # Combine the processed columns with the "Score" column
    paste <(echo -e "Name\ttRNA #\tType\tCodon\tBegin\tEnd\tScore") ${sample}_processed_trnas_temp.tsv ${sample}_score_column.txt > ${sample}_processed_trnas.tsv
    """


}
