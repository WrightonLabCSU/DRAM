process CAT_KEGG_PEP {

    input:
    path( ch_kegg_pep_root_dir )

    output:
    path( "kegg.pep" ), emit: kegg_pep

    script:
    """

    cat ${ch_kegg_pep_root_dir}/*/*pep > kegg.pep

    """
}
