process DBCAN_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file )
    val( top_hit )

    output:
    tuple val( sample ), path ( "${sample}_formatted_hits.out" ), emit: formatted_hits


    output:
    path ( "formatted_dbcan_hits.out" ), emit: formatted_hits

    script:
    """
    python /home/rwoyda/Projects/DRAM2-Nextflow/DRAM2-NF/assets/dbcan_hmm_formatter.py --hits_csv ${hits_file} --db_name ${database_name} --db_handler ${db_handler_file} --output "formatted_dbcan_hits.out"
    
    """
}
