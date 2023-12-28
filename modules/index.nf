process INDEX {

    input:
    path( 'database/*' )
    
    output:
    path( "targetDB.idx" ), emit: index


    script:

    """
    mkdir tmp

    # make query to target db
    mmseqs createindex database/* tmp --threads ${params.threads}

    """
}

