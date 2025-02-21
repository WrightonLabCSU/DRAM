process KOFAM_HMM_FORMATTER {

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"

    tag { sample }

    input:
    tuple val( sample ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    file( ch_kofam_list )

    output:
    tuple val( sample ), path ( "${sample}_formatted_kofam_hits.csv" ), emit: kofam_formatted_hits


    script:
    """
    kofam_hmm_formatter.py --hits_csv ${hits_file} --ch_kofam_ko ${ch_kofam_list} --gene_locs "gene_locs.tsv" --output "${sample}_formatted_kofam_hits.csv"
    
    """
}
