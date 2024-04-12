process COUNT_ANNOTATIONS {

    input:
    file( ch_combined_annotations )
    file( ch_count_annots_script )
    file( ch_distill_sql_script )

    output:
    path "target_id_counts.tsv", emit: target_id_counts
    path "annotations.db", emit: annotations_sqlite3

    script:
    """
    python ${ch_count_annots_script} ${ch_combined_annotations} "target_id_counts.tsv"

    python ${ch_distill_sql_script} --combined_annotations ${ch_combined_annotations} --db_name "annotations.db" 
    """
}
