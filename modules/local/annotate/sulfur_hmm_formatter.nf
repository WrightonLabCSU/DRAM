process SULFUR_HMM_FORMATTER {

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    
    tag { sample }

    input:
    tuple val( sample ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    file( ch_sulfur_formatter )

    output:
    tuple val( sample ), path ( "${sample}_formatted_sulfur_hits.csv" ), emit: sulfur_formatted_hits


    script:
    """
    sulfur_hmm_formatter.py --hits_csv ${hits_file} --gene_locs ${prodigal_locs_tsv} --output "${sample}_formatted_sulfur_hits.csv"
    
    """
}
