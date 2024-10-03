process FORMAT_KEGG_DB {

    input:
    path( ch_kegg_pep, stageAs: "kegg_pep_loc" )
    path( gene_ko_link_loc )
    path( ch_format_kegg_db_script )
    val(kegg_download_date)

    output:
    path( "kegg/*" )

    script:
    """
    if [ ! -f ${gene_ko_link_loc} ]; then
        echo "No Gene KO Link file provided. Running KEGG DB formatting without"
        python ${ch_format_kegg_db_script} --kegg_loc ${ch_kegg_pep} --download_date ${kegg_download_date} --threads ${params.threads} --output_dir kegg
    else
        echo "Gene KO Link file provided. Running KEGG DB formatting with Gene KO Link file"
        python ${ch_format_kegg_db_script} --kegg_loc ${ch_kegg_pep} --gene_ko_link_loc ${gene_ko_link_loc} --download_date ${kegg_download_date} --threads ${params.threads} --output_dir kegg
    fi

    """
}
