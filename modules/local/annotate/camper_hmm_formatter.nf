process CAMPER_HMM_FORMATTER {

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"

    tag { sample }

    input:
    tuple val( sample ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    file( ch_camper_list )

    output:
    tuple val( sample ), path ( "${sample}_formatted_camper_hits.csv" ), emit: camper_formatted_hits


    script:
    """
    camper_hmm_formatter.py --hits_csv ${hits_file} --ch_camper_list ${ch_camper_list} --gene_locs ${prodigal_locs_tsv} --output "${sample}_formatted_camper_hits.csv"
    
    """
}
