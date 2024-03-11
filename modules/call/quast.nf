process QUAST {

    errorStrategy 'finish'

    tag { sample }

    input:
    path ( collected_fasta )
    path ( collected_gff )

    output:
    tuple val( sample ), path( "quast_results/report.tsv" ), emit: quast_tsv
    tuple val( sample ), path( "${sample}_QUAST/icarus.tsv" ), emit: quast_tsv
    tuple val( sample ), path( "${sample}_QUAST/report.html" ), emit: quast_tsv
    tuple val( sample ), path( "${sample}_QUAST/report.pdf" ), emit: quast_tsv

    script:
    """
    quast.py \\
    *.fa \\ 
    -o "quast_results" \\
    --threads ${params.threads}

    mv quast_results/quast_results/* .


    # This will need to be in a loop no to iterate through all .gff files
    # Count number of called genes in ${gff}
    num_genes=\$(grep -c -P "\tgene\t" ${gff})

    # Add a new row to "${sample}_QUAST/${sample}_icarus.tsv" with gene count
    echo "# predicted genes ${num_genes}" >> "${sample}_QUAST/${sample}_icarus.tsv"

    """

}