process COUNT_ANNOTATIONS {
    
    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_biopython:7df21d027f67112e"
    
    input:
    file( ch_combined_annotations )

    output:
    path "target_id_counts.tsv", emit: target_id_counts
    path "annotations.db", emit: annotations_sqlite3

    script:
    """
    count_annotations.py ${ch_combined_annotations} "target_id_counts.tsv"

    distill_sql.py --combined_annotations ${ch_combined_annotations} --db_name "annotations.db" 
    """
}
