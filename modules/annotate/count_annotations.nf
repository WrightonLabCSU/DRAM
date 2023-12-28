process COUNT_ANNOTATIONS {

    input:
    file( combined_annotations )

    output:
    path "target_id_counts.tsv", emit: target_id_counts

    script:
    """
    python /home/rwoyda/Projects/DRAM2-Nextflow/DRAM2-NF/assets/count_annotations.py "${combined_annotations}"
    """
}
