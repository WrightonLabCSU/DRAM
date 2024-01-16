process TRNA_SCAN {

    errorStrategy 'finish'

    tag { sample }

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_processed_trnas.tsv"), emit: trna_scan_out, optional: true

    script:

    """
    tRNAscan-SE \
    -G \
    --thread ${params.threads} \
    -o ${sample}_trna_out.txt \
    ${fasta}

    # Process tRNAscan-SE Output
    if [ -s ${sample}_trna_out.txt ]; then
        trna_frame=\$(awk -F'rsity ' 'NR > 2 { print \$0 }' ${sample}_trna_out.txt | cut -f 1,3-)

        # Remove columns like "Begin.1" or "End.1"
        trna_frame=\$(awk 'BEGIN{OFS=FS="--"}{gsub(/Begin.1|End.1/, "", \$0)}1' <<< "\$trna_frame")

        # Remove duplicate columns
        trna_frame=\$(awk '{a[\$0]++}END{for (i in a) print i}' <<< "\$trna_frame")

        # Insert a new column for the sample name
        trna_frame="${sample} \$trna_frame"

        # Save the processed tRNA frame to a new file
        echo -e "\$trna_frame" > ${sample}_processed_trnas.tsv
    else
        echo "No tRNAs were detected, no trnas.tsv file will be created."
        exit 1
    fi
    """

}
