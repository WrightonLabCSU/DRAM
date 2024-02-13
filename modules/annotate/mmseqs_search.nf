process MMSEQS_SEARCH {

    tag { sample }

    input:
    tuple val( sample ), path( query_database, stageAs: "query_database/" )
    path( mmseqs_database )
    val( bit_score_threshold)
    val( db_name )

    output:
    tuple val( sample ), path("mmseqs_out/${sample}_mmseqs_${db_name}.tsv"), emit: mmseqs_search_out


    script:
    """
    ln -s ${mmseqs_database}/* ${db_name}/

    # Create temporary directory
    mkdir -p mmseqs_out/tmp

    echo "Query Database: ${query_database}"
    echo "Target Database: ${mmseqs_database}"


    # Perform search
    mmseqs search query_database/${sample}.mmsdb ${db_name}/${db_name}.mmsdb mmseqs_out/${sample}_${db_name}.mmsdb mmseqs_out/tmp --threads ${params.threads}

    # Filter to only best hit
    #mmseqs filterdb mmseqs_out/${sample}_${db_name}.mmsdb mmseqs_out/${sample}_${db_name}_tophit.mmsdb --extract-lines 1

    # Filter to only hits with minimum bit score
    #mmseqs filterdb --filter-column 2 --comparison-operator ge --comparison-value ${bit_score_threshold} --threads ${params.threads} mmseqs_out/${sample}_${db_name}_tophit.mmsdb mmseqs_out/${sample}_${db_name}_tophit_minbitscore${bit_score_threshold}.mmsdb

    # Convert results to BLAST outformat 6
    #mmseqs convertalis ${query_database} ${mmseqs_database}  mmseqs_out/${sample}_${db_name}_tophit_minbitscore${bit_score_threshold}.mmsdb mmseqs_out/${sample}_mmseqs_${db_name}.tsv --threads ${params.threads}
    """
}
