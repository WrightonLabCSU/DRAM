process MMSEQS_SEARCH {

    tag { sample }

    input:
    tuple val( sample ), path( query_database, stageAs: "query_database/" )
    path( mmseqs_database )
    val( bit_score_threshold)
    val( db_name )

    output:
    tuple val( sample ), path("mmseqs_out/${sample}_mmseqs_${db_name}.tsv"), emit: mmseqs_search_out
    tuple val( sample ), path("mmseqs_out/${sample}_mmseqs_${db_name}_formatted.tsv"), emit: mmseqs_search_formatted_out

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
    output_path="mmseqs_out/${sample}_mmseqs_${db_name}_formatted.tsv"

    # Use awk to process the file and reorder the columns
    awk -v db_name="${db_name}" 'BEGIN { OFS="\t"; print "query_id", "start_position", "end_position", db_name "_id", db_name "_bitScore" }
    {
        query_id=\$1
        start_position=\$9
        end_position=\$10
        target_id=\$2
        bitScore=\$12
        print query_id, start_position, end_position, target_id, bitScore
    }' "\${input_path}" > "\${output_path}"

    """
}
