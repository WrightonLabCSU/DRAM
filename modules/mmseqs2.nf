process MMSEQS2 {
   // not for use with PFAM
    tag { sample }

    input:
    tuple val( sample ), path( fasta )
    path( 'database/*' )
    path( 'index/*' )

    output:


    // Things we deed defined
    // bit_score_threshold
    // query database - is this the sample
    //  This is, for example, for gene.mmsdb of 
    // target_db -  this the kegg db
    // query target db is the output - needs name
    // query_target_db_top_filt - output alignment file
    // query_target_db_top - output targets to min threshold
    // main output in blast outformat 6 - forward_output_loc


    script:

    """
    mkdir tmp
    

    # make fasta query into database
    # Creates queryDB, queryDB_h, queryDB.index, queryDB_h.index and queryDB.lookup
    mmseqs createdb ${fasta} queryDB

    # make query to target db
    
    mmseqs search queryDB index/* query_target_db tmp --threads ${params.threads}

    """
}


/*
    # filter query to target db to only best hit
    mmseqs \\
    filterdb \\
    ${database} \\ #target
    query_target_db_top \\
    --extract-lines 1 \\

    # filter query to target db to only hits with min threshold
    mmseqs \\
    filterdb \\
    query_target_db_top \\ # output
    query_target_db_top_filt \\
    --filter-column 2 \\
    --comparison-operator e \\
    --comparison-value bit_score_threshold \\
    --threads ${params.threads}

    # convert results to blast outformat 6
    mmseqs \\
    convertalis \\
    ${fasta} \\ #query
    ${database} \\ #target
    query_target_db_top_filt \\ # output alignment file
    forward_output_loc \\
    --threads ${params.threads}
*/