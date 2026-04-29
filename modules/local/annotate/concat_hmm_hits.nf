process CONCAT_HMM_HITS {
    label 'process_tiny'

    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( chunk_csvs, stageAs: "chunk_*.csv" )
    val ( db_name )

    output:
    tuple val( input_fasta ), path( "${input_fasta}___formatted_${db_name}_hits.csv" ), emit: combined_hits

    script:
    """
    out="${input_fasta}___formatted_${db_name}_hits.csv"
    first=1
    for f in chunk_*.csv; do
        if [ "\$first" -eq 1 ]; then
            cat "\$f" > "\$out"
            first=0
        else
            tail -n +2 "\$f" >> "\$out"
        fi
    done
    """
}
