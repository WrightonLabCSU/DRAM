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
    if ! grep -q "^NULL\$" ${db_descriptions}; then
        # Sort the MMseqs output by the column to join on (assuming it's the fourth column, which is `${db_name}_id`)
        cut -d',' -f4 "\${output_path}" | grep -v "${db_name}_id" | sort > "\${output_path}.ids.sorted"
        sort -t$'\t' -k1,1 "${db_descriptions}" > "${db_descriptions}.sorted"

        # Perform the join operation
        # The join field is the first field in db_descriptions and the sorted ids from the MMseqs output
        join -1 1 -2 1 -t$'\t' "${db_descriptions}.sorted" "\${output_path}.ids.sorted" > "\${output_path}.tojoin"

        # Combine the sorted original output with the join results and format it using awk
        awk -v db_name="${db_name}" 'BEGIN { FS=","; OFS="," }
        FNR==NR { a[\$1] = \$0; next }
        {
            split(a[\$1], arr, ",")
            print arr[1], arr[2], arr[3], arr[4], arr[5], \$2, \$3, \$4, \$5, \$6
        }' "\${output_path}" "\${output_path}.tojoin" > "\${output_path}.joined"

        # Move the final joined file to the intended output path
        mv "\${output_path}.joined" "\${output_path}"
    else
        echo "Annotations file contains 'NULL', skipping join operation."
    fi
    
    """
}
