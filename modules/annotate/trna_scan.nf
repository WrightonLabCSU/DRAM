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
    if [ -s ${sample}_trna_out.txt ]; then
        awk -F'\t' 'NR > 2 && !/^-+/{ if (NR == 3) print "Name\ttRNA #\tBegin\tEnd\tType\tCodon\tScore"; if (!/Begin.[0-9]+|End.[0-9]+/) print "${sample}"\t"\$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$6 }' ${sample}_trna_out.txt > ${sample}_processed_trnas.tsv
    else
        echo "No tRNAs were detected, no trnas.tsv file will be created."
        exit 1
    fi
    """

}
