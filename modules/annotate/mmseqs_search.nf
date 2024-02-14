process MMSEQS_SEARCH {

    tag { sample }

    input:
    tuple val( sample ), path( query_database, stageAs: "query_database/" )
    path( mmseqs_database )
    val( bit_score_threshold)
    file( db_descriptions )
    val( db_name )

    output:
    tuple val( sample ), path("mmseqs_out/${sample}_mmseqs_${db_name}.tsv"), emit: mmseqs_search_out
    tuple val( sample ), path("mmseqs_out/${sample}_mmseqs_${db_name}_formatted.csv"), emit: mmseqs_search_formatted_out

    script:
    """
    ln -s ${mmseqs_database}/* ./

    # Create temporary directory
    mkdir -p mmseqs_out/tmp

    echo "Query Database: ${query_database}"
    echo "Target Database: ${mmseqs_database}"


    # Perform search
    mmseqs search query_database/${sample}.mmsdb ${db_name}.mmsdb mmseqs_out/${sample}_${db_name}.mmsdb mmseqs_out/tmp --threads ${params.threads}

    # Filter to only best hit
    mmseqs filterdb mmseqs_out/${sample}_${db_name}.mmsdb mmseqs_out/${sample}_${db_name}_tophit.mmsdb --extract-lines 1

    # Filter to only hits with minimum bit score
    mmseqs filterdb --filter-column 2 --comparison-operator ge --comparison-value ${bit_score_threshold} --threads ${params.threads} mmseqs_out/${sample}_${db_name}_tophit.mmsdb mmseqs_out/${sample}_${db_name}_tophit_minbitscore${bit_score_threshold}.mmsdb

    # Convert results to BLAST outformat 6
    mmseqs convertalis query_database/${sample}.mmsdb ${db_name}.mmsdb  mmseqs_out/${sample}_${db_name}_tophit_minbitscore${bit_score_threshold}.mmsdb mmseqs_out/${sample}_mmseqs_${db_name}.tsv --threads ${params.threads}
    
    # Define the input and output file paths
    input_path="mmseqs_out/${sample}_mmseqs_${db_name}.tsv"
    output_path="mmseqs_out/${sample}_mmseqs_${db_name}_formatted.csv"

    # Use awk to process the file and reorder the columns
    awk -v db_name="${db_name}" 'BEGIN { OFS=","; print "query_id", "start_position", "end_position", db_name "_id", db_name "_bitScore" }
    {
        query_id=\$1
        start_position=\$7
        end_position=\$8
        target_id=\$2
        bitScore=\$12
        print query_id, start_position, end_position, target_id, bitScore
    }' "\${input_path}" > "\${output_path}"

    # Check if the annotations TSV content is not "NULL"
    if ! grep -q "^NULL$" ${db_descriptions}; then
        # Sort the MMseqs output and the additional descriptions file
        sort -k1,1 "\${output_path}" > "\${output_path}.sorted"
        sort -k1,1 "${db_descriptions}" > "${db_descriptions}.sorted"

        # Join the sorted files on the first column, which is the target_id
        # Then use awk to format the output with new header names and concatenate db_name
        join -t, -1 4 -2 1 "\${output_path}.sorted" "${db_descriptions}.sorted" | \
        awk -v db_name="${db_name}" 'BEGIN { OFS=","; print "query_id", "start_position", "end_position", db_name "_id", db_name "_bitScore", db_name "_A_rank", db_name "_B_rank", db_name "_score_type", db_name "_definition", db_name "_ID_for_distillate" }
        {
            # Print the existing fields from MMseqs output
            for (i=1; i<=5; i++) {
                printf "%s,", \$i
            }
            # Print the new fields from db_descriptions file, skipping the first field which is the join field
            for (i=2; i<=NF; i++) {
                printf "%s", \$i (i==NF ? "\n" : ",")
            }
        }' > "\${output_path}.joined"
    fi

    # Move the final joined file to the intended output path
    mv "\${output_path}.joined" "\${output_path}"


    """
}
