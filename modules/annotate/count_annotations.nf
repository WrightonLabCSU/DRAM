process COUNT_ANNOTATIONS {

    input:
    file( combined_annotations )
    file( ch_count_annots_script )

    output:
    path "target_id_counts.tsv", emit: target_id_counts

    script:
    """
    python ${ch_count_annots_script} "${combined_annotations}" "target_id_counts.tsv"
    """
}
