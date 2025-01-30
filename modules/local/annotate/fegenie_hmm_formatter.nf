process FEGENIE_HMM_FORMATTER {

    input:
    tag { sample }

    input:
    tuple val( sample ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )

    output:
    tuple val( sample ), path ( "${sample}_formatted_fegenie_hits.csv" ), emit: fegenie_formatted_hits


    script:
    """
    fegenie_hmm_formatter.py --hits_csv ${hits_file} --gene_locs ${prodigal_locs_tsv} --output "${sample}_formatted_fegenie_hits.csv"
    
    """
}
