process FORMAT_KEGG_DB {

    input:
    path( ch_kegg_pep, stageAs: "kegg_pep_loc" )
    // path( gene_ko_link_loc, stageAs: "genes_ko_link.gz")
    path( ch_format_kegg_db_script )
    val(kegg_download_date)


    script:
    """

    python ${ch_format_kegg_db_script} --kegg_loc ${ch_kegg_pep} --download_date ${kegg_download_date} --threads ${params.threads}

    """
}
