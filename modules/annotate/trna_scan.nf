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
        trna_frame=\$(awk -F'\t' 'NR > 1 && !/^-+/{ print \$0 }' ${sample}_trna_out.txt | cut -f 1,2,3,4,5,6,7,8,9)

        # Remove additional occurrences of "Begin" and "End"
        trna_frame=\$(awk 'BEGIN{OFS=FS="\t"}{gsub(/Begin.[0-9]+|End.[0-9]+/, "", \$0)}1' <<< "\$trna_frame")

        # Remove the "Note" column if it exists
        trna_frame=\$(awk 'BEGIN{OFS=FS="\t"}{NF--; print}' <<< "\$trna_frame")

        # Remove duplicate columns
        trna_frame=\$(awk 'BEGIN{OFS=FS="\t"}{a[\$0]++}END{for (i in a) print i}' <<< "\$trna_frame")

        # Save the processed tRNA frame to a new file, including the header
        echo -e "Name\ttRNA #\tBegin\tEnd\tType\tCodon\tScore" > ${sample}_processed_trnas.tsv
        echo -e "\$trna_frame" >> ${sample}_processed_trnas.tsv
    else
        echo "No tRNAs were detected, no trnas.tsv file will be created."
        exit 1
    fi
    """

}
