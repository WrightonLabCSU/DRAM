process MMSEQS_SEARCH {

    errorStrategy 'finish'

    tag { sample }

    input:
    tuple val( sample ), path( query_database, stageAs: "query_database/" ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    path( mmseqs_database )
    val( bit_score_threshold)
    val( rbh_bit_score_threshold )
    path( db_descriptions, stageAs: "db_descriptions.tsv" )
    val( db_name )
    file( ch_add_db_descriptions )

    output:
    tuple val( sample ), path("mmseqs_out/${sample}_mmseqs_${db_name}_formatted.csv"), emit: mmseqs_search_formatted_out, optional: true
    //tuple val( sample ), path("mmseqs_out/${sample}_mmseqs_rbh_${db_name}.tsv "), emit: mmseqs_search_rbh_formatted_out, optional: true

    script:
    """
    ln -s ${mmseqs_database}/* ./

    # Create temporary directory
    mkdir -p mmseqs_out/tmp

    if [ "${db_name}" != "pfam" ]; then
        # Perform search
        mmseqs search query_database/${sample}.mmsdb ${db_name}.mmsdb mmseqs_out/${sample}_${db_name}.mmsdb mmseqs_out/tmp --threads ${params.threads}

        # Filter to only best hit
        mmseqs filterdb mmseqs_out/${sample}_${db_name}.mmsdb mmseqs_out/${sample}_${db_name}_tophit.mmsdb --extract-lines 1

        # Filter to only hits with minimum bit score
        mmseqs filterdb --filter-column 2 --comparison-operator ge --comparison-value ${bit_score_threshold} --threads ${params.threads} mmseqs_out/${sample}_${db_name}_tophit.mmsdb mmseqs_out/${sample}_${db_name}_tophit_minbitscore${bit_score_threshold}.mmsdb

        # Convert results to BLAST outformat 6
        mmseqs convertalis query_database/${sample}.mmsdb ${db_name}.mmsdb  mmseqs_out/${sample}_${db_name}_tophit_minbitscore${bit_score_threshold}.mmsdb mmseqs_out/${sample}_mmseqs_${db_name}.tsv --threads ${params.threads}

        if [ "${db_name}" == "kegg" ]; then
            # Perform Reverse Best Hit search
            mmseqs search ${db_name}.mmsdb query_database/${sample}.mmsdb  mmseqs_out/${sample}_rbh_${db_name}.mmsdb mmseqs_out/tmp --threads ${params.threads}

            # Filter Reverse Best Hit to only best hit
            mmseqs filterdb mmseqs_out/${sample}_rbh_${db_name}.mmsdb mmseqs_out/${sample}_${db_name}_rbh_tophit.mmsdb --extract-lines 1

            # Filter Reverse Best Hit  to only hits with minimum bit score
            mmseqs filterdb --filter-column 2 --comparison-operator ge --comparison-value ${rbh_bit_score_threshold} --threads ${params.threads} mmseqs_out/${sample}_${db_name}_rbh_tophit.mmsdb mmseqs_out/${sample}_${db_name}_tophit_rbh_minbitscore${bit_score_threshold}.mmsdb

            # Convert Reverse Best Hit  results to BLAST outformat 6
            mmseqs convertalis ${db_name}.mmsdb query_database/${sample}.mmsdb mmseqs_out/${sample}_${db_name}_tophit_rbh_minbitscore${bit_score_threshold}.mmsdb mmseqs_out/${sample}_mmseqs_rbh_${db_name}.tsv --threads ${params.threads}
        
            # Need additional processing for KEGG RBH
        
        fi
    elif [ "${db_name}" == "pfam" ]; then
        # Do profile search:
        mmseqs search query_database/${sample}.mmsdb ${db_name}.mmspro mmseqs_out/${sample}_${db_name}.mmsdb mmseqs_out/tmp --threads ${params.threads}

        # Convert results to BLAST outformat 6
        mmseqs convertalis query_database/${sample}.mmsdb ${db_name}.mmspro mmseqs_out/${sample}_mmseqs_${db_name}.tsv --threads ${params.threads}
    fi

    # Check if the mmseqs_out/${sample}_mmseqs_${db_name}.tsv file is empty
    if [ ! -s "mmseqs_out/${sample}_mmseqs_${db_name}.tsv" ]; then
        echo "The file mmseqs_out/${sample}_mmseqs_${db_name}.tsv is empty. Skipping further processing."
    else
        # Call Python script for further processing
        python ${ch_add_db_descriptions} "${sample}" "${db_name}" "db_descriptions.tsv" "${bit_score_threshold}" "gene_locs.tsv"
    fi

    
    """
}

