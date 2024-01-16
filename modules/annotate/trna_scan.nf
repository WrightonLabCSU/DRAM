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
        trna_frame=\$(awk -F'\t' 'NR > 2 && !/^-+/{ for (i=1; i<=NF; i++) if (!/Begin.[0-9]+|End.[0-9]+/) printf "%s%s", \$i, (i==NF?RS:OFS) }' ${sample}_trna_out.txt)

        # Save the processed tRNA frame to a new file, including the header
        echo -e "Name\ttRNA #\tBegin\tEnd\tType\tCodon\tScore" > ${sample}_processed_trnas.tsv
        echo -e "\$trna_frame" | column -t -s $'\\t' >> ${sample}_processed_trnas.tsv
    else
        echo "No tRNAs were detected, no trnas.tsv file will be created."
        exit 1
    fi
    """

}
