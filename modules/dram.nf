process DRAM {

    input:
    path( clustered_fasta )
    
    output:
    path( "*gene_DB*" ), emit: gene_db
    path( "DRAM_v1.4.4/annotations.tsv" ), emit: annotations
    path( "DRAM_v1.4.4/annotate.log" ), emit: annotations_log
    path( "DRAM_v1.4.4/metabolism_summary.xlsx" ), emit: metab_sum
    path( "DRAM_v1.4.4/distill_log" ), emit: distill_log
    path( "DRAM_v1.4.4/product.html" ), emit: product_html
    path( "DRAM_v1.4.4/product.tsv" ), emit: product


    script:

    """
    source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

    DRAM.py annotate -i ${clustered_fasta} -o DRAM_v1.4.4 --threads ${params.threads}

    DRAM.py distill -i DRAM_v1.4.4/annotations.tsv -o DRAM_v1.4.4/
    
    """
}  