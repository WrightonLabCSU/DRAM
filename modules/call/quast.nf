process QUAST {

    errorStrategy 'finish'

    input:
    path ( collected_fasta )
    path ( collected_gff )

    output:
    path( "quast_results/report.tsv" ), emit: quast_tsv
    path( "${sample}_QUAST/icarus.tsv" )
    path( "${sample}_QUAST/report.html" )
    path( "${sample}_QUAST/report.pdf" )
    path( "collected_quast.tsv" ), emit: quast_collected_out

    script:
    """
    quast.py \\
    *.fa \\ 
    -o "quast_results" \\
    --threads ${params.threads}

    mv quast_results/quast_results/* .


    # Initialize an empty file to collect gene counts
    echo -e "Sample\tPredicted Genes" > gene_counts.tsv

    # Loop through all GFF files and count predicted genes
    for gff in ${collected_gff}; do
        sample_name=\$(basename \$gff .gff)
        num_genes=\$(grep -c -P "\tgene\t" \$gff)
        echo -e "\${sample_name}\t\${num_genes}" >> gene_counts.tsv
    done

    # Add the gene count information to the QUAST report
    cat gene_counts.tsv >> ${sample}_QUAST/report.tsv

    """

}