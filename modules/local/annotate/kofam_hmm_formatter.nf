process KOFAM_HMM_FORMATTER {

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:fc59f737a5e0566a"

    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( hits_file ), path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
    file( ch_kofam_list )

    output:
    tuple val( input_fasta ), path ( "${input_fasta}_formatted_kofam_hits.csv" ), emit: kofam_formatted_hits


    script:
    """
    kofam_hmm_formatter.py --hits_csv ${hits_file} --ch_kofam_ko ${ch_kofam_list} --gene_locs "gene_locs.tsv" --output "${input_fasta}_formatted_kofam_hits.csv"
    
    """
}
